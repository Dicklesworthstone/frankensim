//! Shared gas reflection composes with the original phase and energy owners.
use super::*;
use crate::sublimation::gas::{accommodation::advance_surface_with_accommodation, tests as fixture};

fn controls() -> GasChamberSolve {
    GasChamberSolve { surface: fixture::controls(), maximum_surfaces: 16 }
}
fn law(catalog: &ThermochemicalGasModel, carrier: f64, vapor: f64) -> ThermalAccommodation {
    ThermalAccommodation::new(catalog, &[
        (SpeciesId::new("A").unwrap(), carrier), (fixture::species(), vapor),
    ], [200.0, 900.0]).unwrap()
}
fn surface<'a>(id: u64, phase: &'a dyn SublimationLaw, binding: &'a GasSurfaceBinding,
    solid: SolidSurfaceState, caloric: &'a dyn Fn(f64) -> Result<f64, SublimationError>) -> GasChamberSurface<'a> {
    GasChamberSurface { id, law: phase, binding, solid, temperature_from_enthalpy: caloric }
}

#[test]
fn one_accommodating_chamber_surface_matches_the_original_scalar_endpoint() {
    let catalog = fixture::catalog(-1e7);
    let binding = fixture::binding(-1e7);
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let solid = fixture::solid(400.0, chart);
    let gas = fixture::gas(&catalog, 350.0, 5000.0);
    let accommodation = law(&catalog, 0.5, 0.6);
    let scalar = advance_surface_with_accommodation(&phase, &binding, solid, &gas, &catalog,
        caloric, &accommodation, 0.001, controls().surface, || false).unwrap();
    let maps = BTreeMap::from([(4, accommodation)]);
    let (joint, reflection) = advance_chamber_with_accommodation(
        &[surface(4, &phase, &binding, solid, &caloric)], &maps, &gas, &catalog,
        0.001, controls(), || false,
    ).unwrap();
    assert_eq!(joint.gas_exchange(), scalar.exchange().gas_exchange());
    assert_eq!(joint.surfaces()[0].solid, scalar.exchange().solid());
    assert_eq!(&reflection[&4], scalar.reflection());
    fixture::near(joint.step().energy_residual_j, 0.0, 1e-7);
    fixture::near(joint.step().mass_residual_kg, 0.0, 1e-14);
    let carrier = &joint.step().species[0];
    assert_eq!(carrier.species, SpeciesId::new("A").unwrap());
    assert_eq!(carrier.outgoing_mass_kg, 0.0);
    assert_eq!(carrier.incoming_mass_kg, 0.0);
    assert_eq!(carrier.transferred_mass_kg, 0.0);
    assert!(carrier.energy_into_gas_j > 0.0);
}

#[test]
fn opposing_surfaces_exchange_heat_through_one_unchanged_carrier_inventory() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let cold = fixture::solid(300.0, chart);
    let hot = fixture::solid(500.0, chart);
    let rows = [surface(10, &phase, &binding, cold, &caloric), surface(20, &phase, &binding, hot, &caloric)];
    let maps = BTreeMap::from([(10, law(&catalog, 0.5, 0.0)), (20, law(&catalog, 0.5, 0.0))]);
    let result = advance_chamber_with_accommodation(&rows, &maps, &gas, &catalog, 0.001, controls(), || false).unwrap();
    let reversed = advance_chamber_with_accommodation(&[rows[1], rows[0]], &maps, &gas, &catalog,
        0.001, controls(), || false).unwrap();
    assert_eq!(result, reversed);
    let (candidate, reflection) = result;
    assert_eq!(candidate.gas(), &gas);
    assert_eq!(candidate.step().transferred_mass_kg, 0.0);
    assert_eq!(candidate.step().energy_into_gas_j, 0.0);
    let r = super::super::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
    let expected = 0.5 * 100000.0 / (2.0 * core::f64::consts::PI * r * 400.0).sqrt()
        * 3.0 * r * 100.0 * hot.area_m2 * 0.001;
    fixture::near(reflection[&20].heat_into_gas_j, expected, 1e-12);
    fixture::near(reflection[&10].heat_into_gas_j, -expected, 1e-12);
    fixture::near(candidate.surfaces()[0].surface_temperature_k[1], 300.0 + expected / (cold.mass_kg * chart.1), 1e-9);
    fixture::near(candidate.surfaces()[1].surface_temperature_k[1], 500.0 - expected / (hot.mass_kg * chart.1), 1e-9);
    for row in candidate.surfaces() {
        assert_eq!(row.transferred_mass_kg, 0.0);
        assert_eq!(row.recession_m, 0.0);
    }
    fixture::near(candidate.step().energy_residual_j, 0.0, 1e-7);
}

#[test]
fn joint_phase_and_reflection_use_one_budget_and_keep_specular_default_exact() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let gas = fixture::gas(&catalog, 350.0, 5000.0);
    let rows = [surface(1, &phase, &binding, fixture::solid(500.0, chart), &caloric),
        surface(2, &phase, &binding, fixture::solid(300.0, chart), &caloric)];
    let original = advance_chamber_with_gas(&rows, &gas, &catalog, 0.001, controls(), || false).unwrap();
    let zeros = BTreeMap::from([(1, law(&catalog, 0.0, 0.0)), (2, law(&catalog, 0.0, 0.0))]);
    let (zero, _) = advance_chamber_with_accommodation(&rows, &zeros, &gas, &catalog, 0.001, controls(), || false).unwrap();
    assert_eq!(zero, original);
    let maps = BTreeMap::from([(1, law(&catalog, 0.7, 0.3)), (2, law(&catalog, 0.2, 0.8))]);
    let (joint, reflection) = advance_chamber_with_accommodation(&rows, &maps, &gas, &catalog, 0.001, controls(), || false).unwrap();
    let reflected = reflection.values().map(|row| row.heat_into_gas_j).sum::<f64>();
    fixture::near(joint.step().energy_into_gas_j, original.step().energy_into_gas_j + reflected, 1e-10);
    for (new, old) in joint.surfaces().iter().zip(original.surfaces()) {
        assert_eq!(new.solid.mass_kg, old.solid.mass_kg);
        fixture::near(new.energy_into_gas_j, old.energy_into_gas_j + reflection[&new.id].heat_into_gas_j, 1e-10);
    }
    let a = SpeciesId::new("A").unwrap();
    fixture::near(joint.gas().species_mass_kg(&a).unwrap(), gas.species_mass_kg(&a).unwrap(), 1e-14);
    fixture::near(joint.step().energy_residual_j, 0.0, 1e-7);
    fixture::near(joint.step().mass_residual_kg, 0.0, 1e-14);
    assert!(matches!(advance_chamber_with_accommodation(&rows, &maps, &gas, &catalog,
        2.0, controls(), || false), Err(GasSurfaceError::InventoryLimit { .. })));
}

#[test]
fn missing_surface_laws_late_caloric_failure_and_cancellation_leave_every_input_unchanged() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let solid = fixture::solid(400.0, chart);
    let gas = fixture::gas(&catalog, 350.0, 0.0);
    let saved = gas.clone();
    let rows = [surface(1, &phase, &binding, solid, &caloric), surface(2, &phase, &binding, solid, &caloric)];
    let incomplete = BTreeMap::from([(1, law(&catalog, 0.5, 0.0))]);
    assert!(advance_chamber_with_accommodation(&rows, &incomplete, &gas, &catalog, 0.001, controls(), || false).is_err());
    let extra = BTreeMap::from([(1, law(&catalog, 0.5, 0.0)), (3, law(&catalog, 0.5, 0.0))]);
    assert!(advance_chamber_with_accommodation(&rows, &extra, &gas, &catalog, 0.001, controls(), || false).is_err());
    let maps = BTreeMap::from([(1, law(&catalog, 0.5, 0.0)), (2, law(&catalog, 0.5, 0.0))]);
    let broken = |h| if h == solid.specific_enthalpy_j_kg { caloric(h) } else { Err(SublimationError::Invalid("late caloric domain exit")) };
    let bad = [rows[0], surface(2, &phase, &binding, solid, &broken)];
    assert!(advance_chamber_with_accommodation(&bad, &maps, &gas, &catalog, 0.001, controls(), || false).is_err());
    let mut calls = 0;
    advance_chamber_with_accommodation(&rows, &maps, &gas, &catalog, 0.001, controls(), || { calls += 1; false }).unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(advance_chamber_with_accommodation(&rows, &maps, &gas, &catalog,
            0.001, controls(), || { count += 1; count == stop }).is_err());
        assert_eq!(gas, saved);
        assert_eq!(rows[0].solid, solid);
        assert_eq!(rows[1].solid, solid);
    }
}
