use super::*;
use crate::sublimation::gas::tests as fixture;

fn controls() -> GasChamberAdvanceSolve {
    GasChamberAdvanceSolve {
        step: GasChamberSolve {
            surface: fixture::controls(),
            maximum_surfaces: 8,
        },
        absolute_mass_tolerance_kg: 1e-10,
        absolute_temperature_tolerance_k: 2e-3,
        relative_tolerance: 1e-6,
        minimum_dt_s: 1e-8,
        maximum_dt_s: 1.0,
        maximum_trials: 20000,
    }
}
fn row<'a>(
    id: u64,
    law: &'a dyn fs_thermal::sublimation::SublimationLaw,
    binding: &'a super::super::super::GasSurfaceBinding,
    solid: SolidSurfaceState,
    temperature: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
) -> GasChamberSurface<'a> {
    GasChamberSurface {
        id,
        law,
        binding,
        solid,
        temperature_from_enthalpy: temperature,
    }
}
fn distance(a: &GasChamberCandidate<'_>, b: &GasChamberCandidate<'_>) -> f64 {
    let mut value =
        (a.gas().state().state().temperature - b.gas().state().state().temperature).abs();
    for (a, b) in a.surfaces().iter().zip(b.surfaces()) {
        value += (a.surface_temperature_k[1] - b.surface_temperature_k[1]).abs();
        value += 1e5 * (a.solid.mass_kg - b.solid.mass_kg).abs();
    }
    value
}
fn fixed_reference<'catalog>(
    templates: &[GasChamberSurface<'_>],
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration: f64,
    steps: usize,
    controls: GasChamberSolve,
) -> GasChamberCandidate<'catalog> {
    let mut states: Vec<_> = templates.iter().map(|row| row.solid).collect();
    let mut current = gas.clone();
    let mut endpoint = None;
    for _ in 0..steps {
        let trial = advance_chamber_with_gas(
            &with_states(templates, &states),
            &current,
            catalog,
            duration / steps as f64,
            controls,
            || false,
        )
        .unwrap();
        for (state, row) in states.iter_mut().zip(trial.surfaces()) {
            *state = row.solid;
        }
        current = trial.gas().clone();
        endpoint = Some(trial);
    }
    endpoint.unwrap()
}

#[test]
fn tighter_physical_tolerances_approach_dense_reference_and_close_whole_history() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let law = fixture::phase();
    let temperature = |h| fixture::temperature(h, chart);
    let surfaces = [
        row(
            1,
            &law,
            &binding,
            fixture::solid(420.0, chart),
            &temperature,
        ),
        row(
            2,
            &law,
            &binding,
            fixture::solid(380.0, chart),
            &temperature,
        ),
    ];
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let mut loose_controls = controls();
    loose_controls.absolute_mass_tolerance_kg = 1e-7;
    loose_controls.absolute_temperature_tolerance_k = 0.05;
    loose_controls.relative_tolerance = 1e-4;
    let loose =
        advance_chamber_adaptive(&surfaces, &gas, &catalog, 2.0, loose_controls, || false).unwrap();
    let mut tight_controls = controls();
    tight_controls.absolute_mass_tolerance_kg = 1e-12;
    tight_controls.absolute_temperature_tolerance_k = 1e-4;
    tight_controls.relative_tolerance = 1e-8;
    let tight =
        advance_chamber_adaptive(&surfaces, &gas, &catalog, 2.0, tight_controls, || false).unwrap();
    let reference = fixed_reference(&surfaces, &gas, &catalog, 2.0, 2048, controls().step);
    assert!(distance(tight.endpoint(), &reference) < distance(loose.endpoint(), &reference) * 0.6);
    assert!(tight.audit().accepted_intervals > loose.audit().accepted_intervals);
    assert_eq!(tight.audit().elapsed_s, 2.0);
    assert!(tight.audit().maximum_accepted_error_ratio <= 1.0);
    assert!(tight.audit().surfaces[0].transferred_mass_kg > 0.0);
    assert!(tight.audit().surfaces[1].transferred_mass_kg < 0.0);
    assert!(tight.audit().energy_residual_j.abs() < controls().step.surface.gas.energy_tolerance_j);
    assert!(tight.audit().mass_residual_kg.abs() < controls().step.surface.gas.mass_tolerance_kg);
    for (original, history) in surfaces.iter().zip(&tight.audit().surfaces) {
        let final_mass = tight
            .endpoint()
            .surfaces()
            .iter()
            .find(|row| row.id == history.id)
            .unwrap()
            .solid
            .mass_kg;
        fixture::near(
            original.solid.mass_kg - final_mass,
            history.transferred_mass_kg,
            1e-14,
        );
    }
}

#[test]
fn gross_inventory_limit_refines_without_clipping_or_spending_vapor_twice() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let law = fixture::phase();
    let temperature = |h| fixture::temperature(h, chart);
    let surfaces = [
        row(
            1,
            &law,
            &binding,
            fixture::solid(390.0, chart),
            &temperature,
        ),
        row(
            2,
            &law,
            &binding,
            fixture::solid(380.0, chart),
            &temperature,
        ),
    ];
    let gas = fixture::gas(&catalog, 400.0, 2000.0);
    let mut controls = controls();
    controls.step.surface.maximum_gross_inventory_fraction = 0.01;
    controls.maximum_dt_s = 0.5;
    let result =
        advance_chamber_adaptive(&surfaces, &gas, &catalog, 0.5, controls, || false).unwrap();
    assert!(result.audit().rejected_trials > 0);
    assert!(result.audit().accepted_intervals > 1);
    assert!(
        result
            .audit()
            .surfaces
            .iter()
            .all(|surface| surface.transferred_mass_kg < 0.0)
    );
    assert!(result.gas().species_mass_kg(&fixture::species()).unwrap() > 0.0);
    assert_eq!(result.audit().elapsed_s, 0.5);
    assert!(result.audit().kernel_evaluations <= result.audit().trials * 3);
}

#[test]
fn surface_permutation_preserves_history_and_retry_decisions() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let law = fixture::phase();
    let temperature = |h| fixture::temperature(h, chart);
    let a = row(
        9,
        &law,
        &binding,
        fixture::solid(420.0, chart),
        &temperature,
    );
    let b = row(
        2,
        &law,
        &binding,
        fixture::solid(380.0, chart),
        &temperature,
    );
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let forward =
        advance_chamber_adaptive(&[a, b], &gas, &catalog, 0.75, controls(), || false).unwrap();
    let reverse =
        advance_chamber_adaptive(&[b, a], &gas, &catalog, 0.75, controls(), || false).unwrap();
    assert_eq!(forward, reverse);
    assert_eq!(forward.audit().surfaces[0].id, 2);
    assert_eq!(forward.audit().surfaces[1].id, 9);
}

#[test]
fn trace_vapor_mass_controls_refinement_even_with_loose_temperature_scales() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let law = fixture::phase();
    let temperature = |h| fixture::temperature(h, chart);
    let surfaces = [row(
        1,
        &law,
        &binding,
        fixture::solid(420.0, chart),
        &temperature,
    )];
    let gas = fixture::gas(&catalog, 400.0, 1e-6);
    let mut permissive = controls();
    permissive.absolute_temperature_tolerance_k = 100.0;
    permissive.absolute_mass_tolerance_kg = 1e-3;
    permissive.relative_tolerance = 0.0;
    permissive.maximum_dt_s = 0.1;
    let loose =
        advance_chamber_adaptive(&surfaces, &gas, &catalog, 0.1, permissive, || false).unwrap();
    permissive.absolute_mass_tolerance_kg = 1e-11;
    let tight =
        advance_chamber_adaptive(&surfaces, &gas, &catalog, 0.1, permissive, || false).unwrap();
    assert!(tight.audit().accepted_intervals > loose.audit().accepted_intervals);
    assert!(tight.audit().rejected_trials > loose.audit().rejected_trials);
    assert!(
        tight.gas().species_mass_kg(&fixture::species()).unwrap()
            > gas.species_mass_kg(&fixture::species()).unwrap()
    );
}

#[test]
fn budgets_minimum_step_and_late_cancellation_preserve_all_originals() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let law = fixture::phase();
    let temperature = |h| fixture::temperature(h, chart);
    let surfaces = [row(
        1,
        &law,
        &binding,
        fixture::solid(420.0, chart),
        &temperature,
    )];
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let original_gas = gas.clone();
    let original_solid = surfaces[0].solid;
    let mut limited = controls();
    limited.maximum_trials = 1;
    limited.step.surface.maximum_gross_inventory_fraction = 0.001;
    assert!(matches!(
        advance_chamber_adaptive(&surfaces, &gas, &catalog, 1.0, limited, || false),
        Err(GasChamberAdvanceError::TrialBudget { .. })
    ));
    limited.minimum_dt_s = 0.5;
    assert!(matches!(
        advance_chamber_adaptive(&surfaces, &gas, &catalog, 1.0, limited, || false),
        Err(GasChamberAdvanceError::MinimumStep { .. })
    ));
    let mut calls = 0;
    advance_chamber_adaptive(&surfaces, &gas, &catalog, 0.5, controls(), || {
        calls += 1;
        false
    })
    .unwrap();
    let mut replay_calls = 0;
    let result = advance_chamber_adaptive(&surfaces, &gas, &catalog, 0.5, controls(), || {
        replay_calls += 1;
        replay_calls >= calls - 2
    });
    assert!(matches!(
        result,
        Err(GasChamberAdvanceError::Surface(GasSurfaceError::Surface(
            SublimationError::Cancelled
        )))
    ));
    assert_eq!(gas, original_gas);
    assert_eq!(surfaces[0].solid, original_solid);
}

#[test]
fn zero_duration_and_invalid_controls_have_explicit_outcomes() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let law = fixture::phase();
    let temperature = |h| fixture::temperature(h, chart);
    let surface = row(
        1,
        &law,
        &binding,
        fixture::solid(420.0, chart),
        &temperature,
    );
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let zero =
        advance_chamber_adaptive(&[surface], &gas, &catalog, 0.0, controls(), || false).unwrap();
    assert_eq!(zero.gas(), &gas);
    assert_eq!(zero.endpoint().surfaces()[0].solid, surface.solid);
    assert_eq!(zero.audit().accepted_intervals, 0);
    assert_eq!(zero.audit().energy_into_gas_j, 0.0);
    assert_eq!(zero.audit().surfaces[0].transferred_mass_kg, 0.0);
    let mut bad = controls();
    bad.absolute_mass_tolerance_kg = f64::NAN;
    assert!(matches!(
        advance_chamber_adaptive(&[surface], &gas, &catalog, 1.0, bad, || false),
        Err(GasChamberAdvanceError::Invalid(_))
    ));
    assert!(
        advance_chamber_adaptive(&[surface, surface], &gas, &catalog, 1.0, controls(), || {
            false
        })
        .is_err()
    );
}
