//! Retained loading/failure history under independently refreshed sources.
#[path = "../../adhesive_temperature/common.rs"]
mod original;

use super::{close, controls, directions, driver, potential};
use fs_contact::normal_patch::adhesive::{
    JkrBranch, SourceJkrSphereContact,
    mixed_mode::{
        SourceJkrMixedMode,
        combined::loading::{JkrCombinedLoading, JkrCombinedLoadingError, JkrCombinedLoadingState},
    },
    temperature::{AdhesiveTemperatures, SourceJkrThermalRefresh},
};
use fs_exec::CancelGate;

fn temperatures(a: f64, b: f64, interface: f64) -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: a,
        surface_b_k: b,
        interface_k: interface,
    }
}
fn source(data: &original::Fixture) -> SourceJkrMixedMode {
    let template = super::source();
    SourceJkrMixedMode::new(
        SourceJkrSphereContact::new(
            data.normal.clone(),
            template.normal_source().effective_radius_m(),
            template.normal_source().interaction_range_m(),
            template.normal_source().envelope(),
        )
        .unwrap(),
        template.envelope(),
    )
    .unwrap()
}
fn refresh(
    source: &SourceJkrMixedMode,
    state: JkrCombinedLoadingState,
    data: &original::Fixture,
    temperatures: AdhesiveTemperatures,
) -> (SourceJkrThermalRefresh, SourceJkrMixedMode) {
    let normal = source.normal_source();
    let delta = state.coordinates()[0];
    let witness = if state.contact().is_some() {
        normal.evaluate_attached(delta).unwrap().state()
    } else {
        normal.detached_state(delta).unwrap()
    };
    let refresh = normal
        .refresh_at_held_indentation(witness, data.cards(), temperatures)
        .unwrap();
    let successor = SourceJkrMixedMode::new(refresh.source().clone(), source.envelope()).unwrap();
    (refresh, successor)
}

#[test]
fn heated_joint_loading_recomputes_one_area_boundary_and_all_four_parameter_stores() {
    let data = original::Fixture::with_options(10.0, 0.05, true, 320.0);
    let source = source(&data);
    let gate = CancelGate::new_clock_free();
    let (a, direction) = directions(&source, 4.0, 8.0);
    let delta = -a * a / source.normal_source().effective_radius_m();
    let loading = driver(&source, delta, direction, &gate);
    let previous = loading
        .advance_to_load_factor(loading.initial_state(0.1, &gate).unwrap(), 0.25, &gate)
        .unwrap()
        .state();
    let (normal, successor) = refresh(&source, previous, &data, temperatures(310.0, 305.0, 310.0));
    let continued = loading
        .continue_at_temperature(previous, &normal, &successor, &gate)
        .unwrap();
    let state = continued.state();
    let report = continued.report();
    assert_eq!(continued.previous_state(), previous);
    assert_eq!(state.coordinates(), previous.coordinates());
    assert_eq!(state.load_factor(), previous.load_factor());
    assert_eq!(state.steps(), previous.steps());
    assert!(state.steps() > 0);
    assert_eq!(state.released_energy_j(), previous.released_energy_j());
    assert_eq!(continued.loading().controls(), loading.controls());
    assert_eq!(report.original_model_identity, source.identity());
    assert_eq!(report.refreshed_model_identity, successor.identity());
    assert_eq!(state.model_identity(), successor.identity());
    assert_eq!(report.before, previous.contact());
    assert_eq!(report.after, state.contact());
    assert_eq!(report.previous_limit, loading.limit());
    assert_eq!(report.refreshed_limit, continued.loading().limit());
    assert_ne!(report.previous_limit, report.refreshed_limit);
    assert_ne!(
        previous.contact().unwrap().contact_radius_m,
        state.contact().unwrap().contact_radius_m
    );
    let independent = driver(&successor, delta, direction, &gate);
    assert_eq!(continued.loading().limit(), independent.limit());
    assert_eq!(
        state.contact(),
        independent
            .initial_state(previous.load_factor(), &gate)
            .unwrap()
            .contact()
    );
    assert_eq!(
        report.evaluations,
        2 + independent.preparation_evaluations()
    );
    assert!(report.evaluations <= controls().maximum_evaluations);
    let before = potential(&source, previous.contact().unwrap());
    let after = potential(&successor, state.contact().unwrap());
    let changes = [
        report.normal_elastic_change_j,
        report.tangential_elastic_change_j,
        report.torsional_elastic_change_j,
        report.surface_energy_change_j,
    ];
    for i in 0..4 {
        assert_ne!(changes[i], 0.0);
        close(changes[i], after[i] - before[i], 3e-12);
    }
    close(report.parameter_work_j, changes.into_iter().sum(), 3e-14);
    close(
        report.parameter_work_j,
        after.into_iter().sum::<f64>() - before.into_iter().sum::<f64>(),
        3e-12,
    );
    assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
    assert_ne!(report.parameter_work_j, normal.parameter_work_j());
    // The independent bulk and interface queries retain their selected claims.
    let material = successor.normal_source().normal_state();
    close(
        material.elastic().surface_a_state().young_modulus_pa(),
        9e6,
        2e-14,
    );
    close(
        material.elastic().surface_b_state().young_modulus_pa(),
        19e6,
        2e-14,
    );
    close(material.adhesion_energy_j_per_m2(), 0.035, 2e-14);
    assert_eq!(
        material.interface_state().state_coordinate("T"),
        Some(310.0)
    );
    for (old, new) in source
        .normal_source()
        .normal_state()
        .interface_state()
        .properties()
        .iter()
        .zip(material.interface_state().properties())
    {
        assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
    }
    let (updated, state, _) = continued.into_parts();
    let endpoint = updated.limit().unwrap().load_factor;
    let stopped = updated
        .advance_to_load_factor(state, 2.0 * endpoint, &gate)
        .unwrap();
    assert_eq!(stopped.state().load_factor(), endpoint);
    assert_eq!(stopped.state().steps(), previous.steps() + 1);
    assert_eq!(stopped.report().released_energy_j, 0.0);
    updated.detach_at_fold(stopped.state(), &gate).unwrap();
}

#[test]
fn no_op_preserves_a_cached_fold_and_reversible_temperature_round_trip_returns_history() {
    let data = original::Fixture::with_options(10.0, 0.05, false, 320.0);
    let source = source(&data);
    let gate = CancelGate::new_clock_free();
    let (a, direction) = directions(&source, 4.0, 8.0);
    let delta = -a * a / source.normal_source().effective_radius_m();
    let loading = driver(&source, delta, direction, &gate);
    let fold = loading
        .initial_state(loading.limit().unwrap().load_factor, &gate)
        .unwrap();
    let (same, same_source) = refresh(&source, fold, &data, temperatures(300.0, 300.0, 300.0));
    let unchanged = loading
        .continue_at_temperature(fold, &same, &same_source, &gate)
        .unwrap();
    assert_eq!(unchanged.state(), fold);
    assert_eq!(unchanged.loading().limit(), loading.limit());
    assert_eq!(
        unchanged.loading().preparation_evaluations(),
        loading.preparation_evaluations()
    );
    assert_eq!(unchanged.report().parameter_work_j, 0.0);
    assert_eq!(unchanged.report().storage_change_j, 0.0);
    assert_eq!(unchanged.report().evaluations, 1);
    assert_eq!(
        unchanged
            .loading()
            .detach_at_fold(unchanged.state(), &gate)
            .unwrap(),
        loading.detach_at_fold(fold, &gate).unwrap()
    );
    let previous = loading
        .advance_to_load_factor(loading.initial_state(-0.15, &gate).unwrap(), 0.2, &gate)
        .unwrap()
        .state();
    let (warm, warm_source) = refresh(&source, previous, &data, temperatures(305.0, 295.0, 310.0));
    let warmed = loading
        .continue_at_temperature(previous, &warm, &warm_source, &gate)
        .unwrap();
    let (back, original) = refresh(
        &warm_source,
        warmed.state(),
        &data,
        temperatures(300.0, 300.0, 300.0),
    );
    let restored = warmed
        .loading()
        .continue_at_temperature(warmed.state(), &back, &original, &gate)
        .unwrap();
    // The way back re-resolves the cards with pinned-claim receipts, so the
    // identities differ from the original's single-claim-only ones; the
    // history and every physical value must return exactly.
    original::same_after_refresh(&restored.state(), &previous);
    close(
        restored.report().parameter_work_j,
        -warmed.report().parameter_work_j,
        2e-14,
    );
    original::same_after_refresh(
        &restored
            .loading()
            .advance_to_load_factor(restored.state(), -0.3, &gate)
            .unwrap(),
        &loading
            .advance_to_load_factor(previous, -0.3, &gate)
            .unwrap(),
    );
}

#[test]
fn detached_refresh_survives_loss_of_all_attached_equilibria_without_new_release_or_rebonding() {
    let data = original::Fixture::with_options(10.0, 0.05, false, 320.0);
    let source = source(&data);
    let gate = CancelGate::new_clock_free();
    let (a, direction) = directions(&source, 4.0, 8.0);
    let delta = -2.0 * a * a / source.normal_source().effective_radius_m();
    let loading = driver(&source, delta, direction, &gate);
    let fold = loading
        .initial_state(loading.limit().unwrap().load_factor, &gate)
        .unwrap();
    let previous = loading.detach_at_fold(fold, &gate).unwrap().state();
    assert!(previous.released_energy_j() > 0.0);
    let (hot, successor) = refresh(&source, previous, &data, temperatures(300.0, 300.0, 320.0));
    assert_eq!(hot.before().branch, JkrBranch::Detached);
    assert_eq!(hot.response().branch, JkrBranch::Detached);
    assert!(successor.normal_source().evaluate_attached(delta).is_err());
    let continued = loading
        .continue_at_temperature(previous, &hot, &successor, &gate)
        .unwrap();
    let state = continued.state();
    assert_eq!(state.contact(), None);
    assert_eq!(state.coordinates(), previous.coordinates());
    assert_eq!(state.load_factor(), previous.load_factor());
    assert_eq!(state.steps(), previous.steps());
    assert_eq!(state.released_energy_j(), previous.released_energy_j());
    assert_eq!(continued.loading().limit(), None);
    assert_eq!(continued.loading().preparation_evaluations(), 0);
    assert_eq!(continued.report().previous_limit, loading.limit());
    assert_eq!(continued.report().refreshed_limit, None);
    assert_eq!(continued.report().evaluations, 1);
    assert_eq!(continued.report().parameter_work_j, 0.0);
    assert_eq!(continued.report().storage_change_j, 0.0);
    assert!(continued.loading().initial_state(0.0, &gate).is_err());
    let free = continued
        .loading()
        .advance_to_load_factor(state, -2.0, &gate)
        .unwrap();
    assert_eq!(free.state().contact(), None);
    assert_eq!(free.report().prescribed_work_j, 0.0);
    assert_eq!(free.report().released_energy_j, 0.0);
    assert_eq!(free.report().evaluations, 0);
    assert_eq!(
        free.state().released_energy_j(),
        previous.released_energy_j()
    );
    assert_eq!(free.state().steps(), previous.steps() + 1);
    let (cold, cold_source) = refresh(
        &successor,
        free.state(),
        &data,
        temperatures(300.0, 300.0, 280.0),
    );
    assert!(cold_source.normal_source().evaluate_attached(delta).is_ok());
    let cooled = continued
        .loading()
        .continue_at_temperature(free.state(), &cold, &cold_source, &gate)
        .unwrap();
    assert_eq!(cooled.state().contact(), None);
    assert_eq!(
        cooled.state().released_energy_j(),
        previous.released_energy_j()
    );
    assert_eq!(cooled.report().parameter_work_j, 0.0);
    assert!(cooled.loading().initial_state(0.0, &gate).is_err());
    assert!(matches!(
        cooled.loading().detach_at_fold(cooled.state(), &gate),
        Err(JkrCombinedLoadingError::NotAtSeparatingFold)
    ));
}

#[test]
fn heated_combined_instability_and_foreign_source_controls_cannot_clamp_or_replace_history() {
    let data = original::Fixture::with_options(10.0, 0.05, false, 320.0);
    let source = source(&data);
    let gate = CancelGate::new_clock_free();
    let (a, direction) = directions(&source, 4.0, 8.0);
    let delta = -a * a / source.normal_source().effective_radius_m();
    let loading = driver(&source, delta, direction, &gate);
    let previous = loading.initial_state(0.85, &gate).unwrap();
    let (hot, successor) = refresh(&source, previous, &data, temperatures(300.0, 300.0, 320.0));
    assert_eq!(hot.response().branch, JkrBranch::Attached);
    assert!(
        driver(&successor, delta, direction, &gate)
            .limit()
            .unwrap()
            .load_factor
            < previous.load_factor()
    );
    assert!(matches!(
        loading.continue_at_temperature(previous, &hot, &successor, &gate),
        Err(JkrCombinedLoadingError::BeyondFirstBoundary)
    ));
    let previous = loading.initial_state(0.2, &gate).unwrap();
    let (warm, successor) = refresh(&source, previous, &data, temperatures(305.0, 295.0, 310.0));
    assert!(matches!(
        loading.continue_at_temperature(previous, &warm, &source, &gate),
        Err(JkrCombinedLoadingError::Invalid(_))
    ));
    let mut envelope = source.envelope();
    envelope.maximum_tangential_displacement_to_radius *= 1.1;
    let changed_envelope = SourceJkrMixedMode::new(warm.source().clone(), envelope).unwrap();
    assert!(matches!(
        loading.continue_at_temperature(previous, &warm, &changed_envelope, &gate),
        Err(JkrCombinedLoadingError::Invalid(_))
    ));
    let foreign = driver(
        &source,
        delta,
        [-direction[0], direction[1], direction[2]],
        &gate,
    );
    assert!(matches!(
        loading.continue_at_temperature(
            foreign.initial_state(0.2, &gate).unwrap(),
            &warm,
            &successor,
            &gate
        ),
        Err(JkrCombinedLoadingError::StateMismatch)
    ));
    let normal = source.normal_source();
    let wrong_approach = normal
        .refresh_at_held_indentation(
            normal.evaluate_attached(0.9 * delta).unwrap().state(),
            data.cards(),
            temperatures(305.0, 295.0, 310.0),
        )
        .unwrap();
    assert!(matches!(
        loading.continue_at_temperature(previous, &wrong_approach, &successor, &gate),
        Err(JkrCombinedLoadingError::Invalid(_))
    ));
    let wrong_branch = normal
        .refresh_at_held_indentation(
            normal.detached_state(delta).unwrap(),
            data.cards(),
            temperatures(305.0, 295.0, 310.0),
        )
        .unwrap();
    assert!(matches!(
        loading.continue_at_temperature(previous, &wrong_branch, &successor, &gate),
        Err(JkrCombinedLoadingError::Invalid(_))
    ));
    let held_load = normal
        .refresh_at_held_load(
            normal.evaluate_attached(0.0).unwrap().normal_force_n,
            data.cards(),
            temperatures(300.0, 300.0, 300.0),
        )
        .unwrap();
    let held_source =
        SourceJkrMixedMode::new(held_load.source().clone(), source.envelope()).unwrap();
    assert!(matches!(
        loading.continue_at_temperature(previous, &held_load, &held_source, &gate),
        Err(JkrCombinedLoadingError::Invalid(_))
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        loading.continue_at_temperature(previous, &warm, &successor, &cancelled),
        Err(JkrCombinedLoadingError::Cancelled)
    ));
    assert_eq!(
        loading
            .initial_state(previous.load_factor(), &gate)
            .unwrap(),
        previous
    );
    loading
        .continue_at_temperature(previous, &warm, &successor, &gate)
        .unwrap();
}

#[test]
fn successor_preparation_and_actual_endpoint_share_one_original_evaluation_budget() {
    let data = original::Fixture::with_options(10.0, 0.05, false, 320.0);
    let source = source(&data);
    let gate = CancelGate::new_clock_free();
    let (_, direction) = directions(&source, 4.0, 0.0);
    for maximum in [2, 3] {
        let loading = JkrCombinedLoading::new(
            &source, 0.0, [0.0; 2], direction[2],
            fs_contact::normal_patch::adhesive::mixed_mode::combined::loading::JkrCombinedLoadingControls {
                maximum_evaluations: maximum,
                ..controls()
            },
            &gate,
        ).unwrap();
        assert_eq!(loading.preparation_evaluations(), 1);
        assert_eq!(loading.limit(), None);
        let previous = loading.initial_state(0.2, &gate).unwrap();
        let (warm, successor) =
            refresh(&source, previous, &data, temperatures(305.0, 295.0, 310.0));
        let candidate = loading.continue_at_temperature(previous, &warm, &successor, &gate);
        if maximum == 2 {
            assert!(matches!(candidate, Err(JkrCombinedLoadingError::Budget)));
        } else {
            let candidate = candidate.unwrap();
            assert_eq!(candidate.report().evaluations, 3);
            assert_eq!(candidate.loading().preparation_evaluations(), 1);
            assert_eq!(candidate.state().steps(), previous.steps());
            assert_eq!(candidate.state().coordinates(), previous.coordinates());
        }
        assert_eq!(loading.initial_state(0.2, &gate).unwrap(), previous);
    }
}
