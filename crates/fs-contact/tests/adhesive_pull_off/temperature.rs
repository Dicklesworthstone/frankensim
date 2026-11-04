//! G1/G3/G4: held apparatus stages with original thermal adhesion and shear.
#[path = "../adhesive_temperature/common.rs"]
mod source_fixture;

use core::f64::consts::PI;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::mixed_mode::loading::JkrShearBoundaryKind;
use fs_contact::normal_patch::adhesive::mixed_mode::pull_off::temperature::*;
use fs_contact::normal_patch::adhesive::mixed_mode::pull_off::*;
use fs_contact::normal_patch::adhesive::mixed_mode::{JkrMixedModeEnvelope, SourceJkrMixedMode};
use fs_contact::normal_patch::adhesive::temperature::AdhesiveTemperatures;
use fs_exec::CancelGate;
use fs_material::state_point::MaterialPropertySelection;

fn temperatures(a: f64, b: f64, interface: f64) -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: a,
        surface_b_k: b,
        interface_k: interface,
    }
}
fn fixture(data: &source_fixture::Fixture) -> JkrPullOffFixture {
    fixture_with_shear(data, [0.3, -0.2])
}
fn fixture_with_shear(data: &source_fixture::Fixture, shear_stage: [f64; 2]) -> JkrPullOffFixture {
    let normal = data.jkr();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let a = (PI * normal.normal_state().adhesion_energy_j_per_m2() / (8.0 * e)).cbrt();
    let depth = a * a;
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    JkrPullOffFixture::new(
        SourceJkrMixedMode::new(
            normal,
            JkrMixedModeEnvelope {
                maximum_tangential_displacement_to_radius: 0.1,
                maximum_mean_shear_stress_to_yield: 0.1,
            },
        )
        .unwrap(),
        shear,
        4.0 * e * a,
        4.0 * e * a,
        shear_stage.map(|x| x * depth),
        JkrPullOffLimits {
            maximum_contact_radius_m: 4.0 * a,
            maximum_absolute_normal_stage_m: 1000.0 * depth,
            maximum_iterations: 96,
            maximum_evaluations: 256,
        },
        &CancelGate::new(),
    )
    .unwrap()
}
fn interior(law: &JkrPullOffFixture) -> JkrPullOffResponse {
    let stage = 0.5
        * (law.initial_response().state().normal_stage_m()
            + law.first_boundary().response.state().normal_stage_m());
    law.evaluate_attached(stage, &CancelGate::new()).unwrap()
}
fn close(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 2e-12 * scale, "{a:e} != {b:e}");
}

#[test]
fn unequal_temperatures_move_the_actual_mixed_mode_tip_and_both_spring_stores() {
    // Hot-only alternatives make an unpinned re-resolution ambiguous.
    let data = source_fixture::Fixture::with_options(1.0, 1.0, true, 320.0);
    let law = fixture(&data);
    let saved = law.clone();
    let before = interior(&law);
    let update = law
        .refresh_at_held_stages(
            before.state(),
            data.cards(),
            temperatures(310.0, 290.0, 315.0),
            &CancelGate::new(),
        )
        .unwrap();
    let after = update.response();
    let next = update.source();
    let normal = next.contact().normal_source().normal_state();
    assert_eq!(update.before(), &before);
    assert_eq!(update.original_fixture_identity(), law.identity());
    assert_eq!(
        after.state().normal_stage_m(),
        before.state().normal_stage_m()
    );
    assert_eq!(after.state().branch(), JkrPullOffBranch::Attached);
    assert_eq!(
        next.normal_stiffness_n_per_m(),
        law.normal_stiffness_n_per_m()
    );
    assert_eq!(
        next.tangential_stiffness_n_per_m(),
        law.tangential_stiffness_n_per_m()
    );
    assert_eq!(next.tangential_stage_m(), law.tangential_stage_m());
    assert_ne!(after.tip_coordinates_m, before.tip_coordinates_m);
    assert_ne!(
        after.contact.unwrap().contact_area_m2,
        before.contact.unwrap().contact_area_m2
    );
    assert_ne!(after.normal_spring_energy_j, before.normal_spring_energy_j);
    assert_ne!(
        after.tangential_spring_energy_j,
        before.tangential_spring_energy_j
    );
    assert_ne!(after.state(), next.initial_response().state());
    close(
        normal.elastic().surface_a_state().young_modulus_pa(),
        9e5,
        1e6,
    );
    close(
        normal.elastic().surface_b_state().young_modulus_pa(),
        2.2e6,
        2e6,
    );
    close(normal.adhesion_energy_j_per_m2(), 0.55, 1.0);
    close(
        next.shear_source().law().static_shear_strength_pa(),
        2550.0,
        3000.0,
    );
    assert_eq!(normal.interface_state().state_coordinate("T"), Some(315.0));
    assert_eq!(normal.interface_state().state_coordinate("p"), Some(2e5));
    assert_eq!(
        next.shear_source().normal_state_identity(),
        normal.identity()
    );
    assert_eq!(
        after.contact.unwrap().normal_state_identity(),
        normal.identity()
    );
    for (old, new) in law
        .shear_source()
        .resolved()
        .properties()
        .iter()
        .zip(next.shear_source().resolved().properties())
    {
        assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
    }
    // External stages do not move; internal tip motion exchanges energy with
    // both springs, so contact-potential change alone is the wrong supply.
    let contact_change =
        after.contact.unwrap().reversible_energy_j - before.contact.unwrap().reversible_energy_j;
    let spring_change = after.normal_spring_energy_j - before.normal_spring_energy_j
        + after.tangential_spring_energy_j
        - before.tangential_spring_energy_j;
    close(
        update.parameter_work_j(),
        contact_change + spring_change,
        0.1,
    );
    assert_eq!(update.mechanical_work_j(), 0.0);
    assert_eq!(update.potential_change_j(), update.parameter_work_j());
    assert_eq!(after.state().accumulated_released_energy_j(), 0.0);
    assert_eq!(law, saved);
    let step = next
        .advance(
            after.state(),
            0.95 * after.state().normal_stage_m(),
            &CancelGate::new(),
        )
        .unwrap();
    assert!(step.event.is_none());
    assert_eq!(step.released_energy_j, 0.0);
    assert_eq!(
        next.advance(
            before.state(),
            after.state().normal_stage_m(),
            &CancelGate::new()
        ),
        Err(JkrPullOffError::StateMismatch)
    );
}

#[test]
fn no_op_is_exact_and_a_supported_thermal_cycle_restores_equilibrium_and_work() {
    let data = source_fixture::Fixture::new();
    let law = fixture(&data);
    let before = interior(&law);
    let same = law
        .refresh_at_held_stages(
            before.state(),
            data.cards(),
            temperatures(300.0, 300.0, 300.0),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(same.source(), &law);
    assert_eq!(same.before(), &before);
    assert_eq!(same.response(), &before);
    assert_eq!(same.parameter_work_j(), 0.0);
    let hot = law
        .refresh_at_held_stages(
            before.state(),
            data.cards(),
            temperatures(305.0, 295.0, 310.0),
            &CancelGate::new(),
        )
        .unwrap();
    let cold = hot
        .source()
        .refresh_at_held_stages(
            hot.response().state(),
            data.cards(),
            temperatures(300.0, 300.0, 300.0),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(cold.source(), &law);
    assert_eq!(cold.response(), &before);
    assert_eq!(hot.parameter_work_j() + cold.parameter_work_j(), 0.0);
    assert_eq!(cold.clone().into_source(), law);
}

#[test]
fn detached_refresh_preserves_released_energy_and_never_rebonds_or_reemits_it() {
    let data = source_fixture::Fixture::new();
    let law = fixture(&data);
    let boundary = law.first_boundary();
    assert!(boundary.relaxed_separation_admissible);
    let target = 2.0 * boundary.response.state().normal_stage_m();
    let separation = law
        .advance(law.initial_response().state(), target, &CancelGate::new())
        .unwrap();
    assert_eq!(
        separation.response.state().branch(),
        JkrPullOffBranch::Detached
    );
    assert!(separation.released_energy_j > 0.0);
    let update = law
        .refresh_at_held_stages(
            separation.response.state(),
            data.cards(),
            temperatures(310.0, 290.0, 315.0),
            &CancelGate::new(),
        )
        .unwrap();
    let after = update.response();
    assert_eq!(update.before(), &separation.response);
    assert_eq!(after.state().branch(), JkrPullOffBranch::Detached);
    assert_eq!(after.state().normal_stage_m(), target);
    assert_eq!(
        after.state().accumulated_released_energy_j(),
        separation.released_energy_j
    );
    assert_eq!(after.contact, None);
    assert_eq!(after.stage_force_n, [0.0; 3]);
    assert_eq!(after.total_stored_energy_j, 0.0);
    assert_eq!(update.parameter_work_j(), 0.0);
    let further = update
        .source()
        .advance(after.state(), 1.1 * target, &CancelGate::new())
        .unwrap();
    assert_eq!(further.released_energy_j, 0.0);
    assert_eq!(further.event, None);
    assert_eq!(
        further.response.state().accumulated_released_energy_j(),
        separation.released_energy_j
    );
}

#[test]
fn lost_attachment_foreign_cards_expired_sources_and_cancel_are_atomic() {
    let data = source_fixture::Fixture::new();
    let law = fixture(&data);
    let saved = law.clone();
    let before = interior(&law);
    let fold = law.first_boundary().response;
    let near = law
        .evaluate_attached(0.98 * fold.state().normal_stage_m(), &CancelGate::new())
        .unwrap();
    assert!(matches!(
        law.refresh_at_held_stages(
            near.state(),
            data.cards(),
            temperatures(300.0, 300.0, 320.0),
            &CancelGate::new(),
        ),
        Err(JkrPullOffTemperatureError::Fixture(
            JkrPullOffError::BeyondFirstBoundary
        ))
    ));
    // A stronger held shear stage loses source-admitted stick before the
    // radius fold. Thermal re-equilibration must not bypass that first event.
    let sheared = fixture_with_shear(&data, [3.0, 0.0]);
    let hot = sheared
        .refresh_at_held_stages(
            interior(&sheared).state(),
            data.cards(),
            temperatures(300.0, 300.0, 305.0),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(
        hot.source().first_boundary().kind,
        JkrShearBoundaryKind::GrossSlip
    );
    let held_stage = 1.02
        * hot
            .source()
            .first_boundary()
            .response
            .state()
            .normal_stage_m();
    assert!(held_stage > sheared.first_boundary().response.state().normal_stage_m());
    let before_slip = sheared
        .evaluate_attached(held_stage, &CancelGate::new())
        .unwrap();
    assert!(matches!(
        sheared.refresh_at_held_stages(
            before_slip.state(),
            data.cards(),
            temperatures(300.0, 300.0, 305.0),
            &CancelGate::new(),
        ),
        Err(JkrPullOffTemperatureError::Fixture(
            JkrPullOffError::BeyondFirstBoundary
        ))
    ));
    assert!(
        law.refresh_at_held_stages(
            before.state(),
            [&data.cards[1], &data.cards[0]],
            temperatures(300.0, 300.0, 300.0),
            &CancelGate::new(),
        )
        .is_err()
    );
    let expired_data = source_fixture::Fixture::with_options(1.0, 1.0, true, 310.0);
    let expired = fixture(&expired_data);
    assert!(matches!(
        expired.refresh_at_held_stages(
            interior(&expired).state(),
            expired_data.cards(),
            temperatures(300.0, 300.0, 315.0),
            &CancelGate::new(),
        ),
        Err(JkrPullOffTemperatureError::Source(_))
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        law.refresh_at_held_stages(
            before.state(),
            data.cards(),
            temperatures(300.0, 300.0, 300.0),
            &cancelled
        ),
        Err(JkrPullOffTemperatureError::Fixture(
            JkrPullOffError::Cancelled
        ))
    );
    assert_eq!(law, saved);
    assert_eq!(
        law.evaluate_attached(before.state().normal_stage_m(), &CancelGate::new())
            .unwrap(),
        before
    );
    assert!(
        law.refresh_at_held_stages(
            before.state(),
            data.cards(),
            temperatures(310.0, 290.0, 315.0),
            &CancelGate::new(),
        )
        .is_ok()
    );
}
