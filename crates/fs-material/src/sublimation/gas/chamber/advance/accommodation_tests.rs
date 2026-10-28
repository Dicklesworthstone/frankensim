//! Numerical integration tests with manufactured, calorically compatible phases.
//! Coefficients are explicit test data, not measured material accommodation.
use super::*;
use crate::sublimation::gas::{GasSurfaceBinding, tests as fixture};
use fs_thermal::sublimation::SublimationLaw;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn controls() -> GasChamberAdvanceSolve {
    GasChamberAdvanceSolve {
        step: GasChamberSolve { surface: fixture::controls(), maximum_surfaces: 8 },
        absolute_mass_tolerance_kg: 1e-11,
        absolute_temperature_tolerance_k: 1e-5,
        relative_tolerance: 0.0,
        minimum_dt_s: 1e-8,
        maximum_dt_s: 0.1,
        maximum_trials: 10000,
    }
}
fn accommodation(catalog: &ThermochemicalGasModel, alpha: f64) -> ThermalAccommodation {
    ThermalAccommodation::new(catalog, &[
        (SpeciesId::new("A").unwrap(), alpha), (fixture::species(), alpha),
    ], [200.0, 900.0]).unwrap()
}
fn row<'a>(id: u64, phase: &'a dyn SublimationLaw, binding: &'a GasSurfaceBinding,
    solid: SolidSurfaceState, caloric: &'a dyn Fn(f64) -> Result<f64, SublimationError>) -> GasChamberSurface<'a> {
    GasChamberSurface { id, law: phase, binding, solid, temperature_from_enthalpy: caloric }
}

#[test]
fn pure_reflection_refines_against_closed_form_without_mass_error_to_drive_it() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let solid = fixture::solid(400.0, chart);
    let gas = fixture::gas(&catalog, 300.0, 0.0);
    let rows = [row(7, &phase, &binding, solid, &caloric)];
    let laws = BTreeMap::from([(7, accommodation(&catalog, 0.5))]);
    let mut loose_controls = controls();
    loose_controls.absolute_temperature_tolerance_k = 0.05;
    let loose = advance_chamber_adaptive_with_accommodation(&rows, &laws, &gas,
        &catalog, 0.1, loose_controls, || false).unwrap();
    let tight = advance_chamber_adaptive_with_accommodation(&rows, &laws, &gas,
        &catalog, 0.1, controls(), || false).unwrap();

    // Fixed masses and ideal-diatomic carrier: q=K sqrt(Tg) (Ts-Tg).
    // Conservation reduces this to dTg/dt = K(1+Cg/Cs)/Cg sqrt(Tg)(Teq-Tg).
    // Integrating its separable scalar equation supplies an independent oracle.
    let r = R / 0.02897;
    let cs = solid.mass_kg * chart.1;
    let cg = gas.mass_kg() * 2.5 * r;
    let equilibrium = (cs * 400.0 + cg * 300.0) / (cs + cg);
    let k = 0.5 * solid.area_m2 * gas.mass_kg() / 0.01
        * (r / (2.0 * core::f64::consts::PI)).sqrt() * 3.0 * r;
    let x = (300.0 / equilibrium).sqrt();
    let z = 0.5 * fs_math::det::ln((1.0 + x) / (1.0 - x))
        + k * (1.0 + cg / cs) * equilibrium.sqrt() * 0.1 / (2.0 * cg);
    let v = fs_math::det::tanh(z);
    let exact = equilibrium * v * v;
    let loose_error = (loose.gas().state().state().temperature - exact).abs();
    let tight_error = (tight.gas().state().state().temperature - exact).abs();
    assert!(tight_error < 0.2 * loose_error, "{tight_error} / {loose_error}");
    assert!(tight.audit().accepted_intervals > loose.audit().accepted_intervals);
    assert!(tight.audit().rejected_trials > 0);
    assert!(tight.audit().maximum_accepted_error_ratio <= 1.0);
    assert_eq!(tight.audit().elapsed_s, 0.1);
    assert_eq!(tight.gas().mass_kg(), gas.mass_kg());
    assert_eq!(tight.endpoint().surfaces()[0].solid.mass_kg, solid.mass_kg);
    assert_eq!(tight.audit().surfaces[0].transferred_mass_kg, 0.0);
    assert_eq!(tight.audit().surfaces[0].recession_m, 0.0);
    assert!(tight.audit().kernel_evaluations <= 3 * tight.audit().trials);
    // This fails if coarse or rejected probe heat leaks into accepted history.
    let heat = cg * (tight.gas().state().state().temperature - 300.0);
    fixture::near(tight.audit().reflected_heat_into_gas_j, heat, 1e-7);
    fixture::near(tight.audit().energy_into_gas_j, heat, 1e-7);
    fixture::near(tight.audit().surfaces[0].reflected_heat_into_gas_j, heat, 1e-7);
    fixture::near(tight.audit().energy_residual_j, 0.0, 1e-7);
    fixture::near(tight.audit().gas_stream_energy_residual_j, 0.0, 1e-7);
    fixture::near(tight.audit().solid_stream_energy_residual_j, 0.0, 1e-7);
}

#[test]
fn zero_coefficients_preserve_specular_states_histories_and_retry_decisions() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let rows = [row(2, &phase, &binding, fixture::solid(420.0, chart), &caloric),
        row(1, &phase, &binding, fixture::solid(380.0, chart), &caloric)];
    let laws = BTreeMap::from([(1, accommodation(&catalog, 0.0)), (2, accommodation(&catalog, 0.0))]);
    let original = advance_chamber_adaptive(&rows, &gas, &catalog, 0.1, controls(), || false).unwrap();
    let extended = advance_chamber_adaptive_with_accommodation(&rows, &laws, &gas,
        &catalog, 0.1, controls(), || false).unwrap();
    assert_eq!(original, extended);
    assert_eq!(extended.audit().reflected_heat_into_gas_j, 0.0);
}

#[test]
fn heterogeneous_coefficients_keep_canonical_order_and_actual_energy_conservation() {
    let catalog = fixture::catalog(-1e7);
    let binding = fixture::binding(-1e7);
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let gas = fixture::gas(&catalog, 350.0, 1000.0);
    let rows = [row(20, &phase, &binding, fixture::solid(420.0, chart), &caloric),
        row(10, &phase, &binding, fixture::solid(300.0, chart), &caloric)];
    let laws = BTreeMap::from([(10, accommodation(&catalog, 0.2)), (20, accommodation(&catalog, 0.7))]);
    let controls = GasChamberAdvanceSolve { maximum_dt_s: 0.001, ..controls() };
    let forward = advance_chamber_adaptive_with_accommodation(&rows, &laws, &gas,
        &catalog, 0.01, controls, || false).unwrap();
    let reverse = advance_chamber_adaptive_with_accommodation(&[rows[1], rows[0]], &laws,
        &gas, &catalog, 0.01, controls, || false).unwrap();
    assert_eq!(forward, reverse);
    assert_eq!(forward.audit().surfaces[0].id, 10);
    assert!(forward.audit().surfaces[0].reflected_heat_into_gas_j < 0.0);
    assert!(forward.audit().surfaces[1].reflected_heat_into_gas_j > 0.0);
    fixture::near(forward.audit().reflected_heat_into_gas_j,
        forward.audit().surfaces.iter().map(|s| s.reflected_heat_into_gas_j).sum(), 1e-10);
    fixture::near(forward.audit().energy_residual_j, 0.0, 1e-7);
    fixture::near(forward.audit().mass_residual_kg, 0.0, 1e-14);
    assert!(forward.audit().maximum_species_mass_residual_kg <= 1e-14);
    fixture::near(forward.gas().species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap(),
        gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap(), 1e-14);
}

#[test]
fn gross_inventory_retries_preserve_reflection_accounting() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 2000.0);
    let rows = [row(1, &phase, &binding, fixture::solid(380.0, chart), &caloric),
        row(2, &phase, &binding, fixture::solid(390.0, chart), &caloric)];
    let laws = BTreeMap::from([(1, accommodation(&catalog, 0.02)), (2, accommodation(&catalog, 0.02))]);
    let mut controls = controls();
    controls.step.surface.maximum_gross_inventory_fraction = 0.001;
    let accepted = advance_chamber_adaptive_with_accommodation(&rows, &laws, &gas,
        &catalog, 0.02, controls, || false).unwrap();
    assert!(accepted.audit().rejected_trials > 0);
    assert!(accepted.audit().reflected_heat_into_gas_j < 0.0);
    assert!(accepted.gas().species_mass_kg(&fixture::species()).unwrap() > 0.0);
    fixture::near(accepted.audit().gas_stream_energy_residual_j, 0.0, 1e-7);
    fixture::near(accepted.audit().solid_stream_energy_residual_j, 0.0, 1e-7);
}

#[test]
fn full_roster_is_required_even_at_zero_duration() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let gas = fixture::gas(&catalog, 300.0, 0.0);
    let rows = [row(1, &phase, &binding, fixture::solid(400.0, chart), &caloric)];
    let laws = BTreeMap::from([(1, accommodation(&catalog, 0.5))]);
    let no_time = advance_chamber_adaptive_with_accommodation(&rows, &laws, &gas,
        &catalog, 0.0, controls(), || false).unwrap();
    assert_eq!(no_time.gas(), &gas);
    assert_eq!(no_time.endpoint().surfaces()[0].solid, rows[0].solid);
    assert_eq!(no_time.audit().reflected_heat_into_gas_j, 0.0);
    assert_eq!(no_time.audit().accepted_intervals, 0);
    assert_eq!(no_time.audit().kernel_evaluations, 1);
    for bad in [BTreeMap::new(), BTreeMap::from([(2, accommodation(&catalog, 0.5))]),
        BTreeMap::from([(1, accommodation(&catalog, 0.5)), (2, accommodation(&catalog, 0.5))])] {
        assert!(advance_chamber_adaptive_with_accommodation(&rows, &bad, &gas,
            &catalog, 0.0, controls(), || false).is_err());
    }
}

#[test]
fn heat_driven_work_limits_and_late_cancellation_publish_no_history() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let gas = fixture::gas(&catalog, 300.0, 0.0);
    let saved = gas.clone();
    let solid = fixture::solid(400.0, chart);
    let rows = [row(1, &phase, &binding, solid, &caloric)];
    let laws = BTreeMap::from([(1, accommodation(&catalog, 0.5))]);
    let limited = GasChamberAdvanceSolve { maximum_trials: 1, ..controls() };
    assert!(matches!(advance_chamber_adaptive_with_accommodation(&rows, &laws, &gas,
        &catalog, 0.1, limited, || false), Err(GasChamberAdvanceError::TrialBudget { .. })));
    let minimum = GasChamberAdvanceSolve { minimum_dt_s: 0.05, ..controls() };
    assert!(matches!(advance_chamber_adaptive_with_accommodation(&rows, &laws, &gas,
        &catalog, 0.1, minimum, || false), Err(GasChamberAdvanceError::MinimumStep { .. })));
    let mut calls = 0;
    advance_chamber_adaptive_with_accommodation(&rows, &laws, &gas,
        &catalog, 0.1, controls(), || { calls += 1; false }).unwrap();
    let mut count = 0;
    assert!(advance_chamber_adaptive_with_accommodation(&rows, &laws, &gas,
        &catalog, 0.1, controls(), || { count += 1; count == calls }).is_err());
    assert_eq!(gas, saved);
    assert_eq!(rows[0].solid, solid);
}
