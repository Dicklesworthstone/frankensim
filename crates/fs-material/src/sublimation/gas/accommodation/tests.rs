//! Manufactured kinetic/caloric limits, not measured accommodation coefficients.
use super::*;
use crate::sublimation::gas::{advance_surface_with_gas, tests as fixture};
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

fn accommodation(catalog: &ThermochemicalGasModel, a: f64, b: f64) -> ThermalAccommodation {
    ThermalAccommodation::new(catalog, &[
        (SpeciesId::new("B").unwrap(), b), (SpeciesId::new("A").unwrap(), a),
    ], [200.0, 900.0]).unwrap()
}

#[test]
fn carrier_heats_without_mass_transfer_and_uses_the_original_finite_caloric_stores() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let chart = fixture::chart(&catalog, &binding);
    let solid = fixture::solid(400.0, chart);
    let gas = fixture::gas(&catalog, 300.0, 0.0);
    let law = accommodation(&catalog, 0.5, 0.0);
    let dt = 0.001;
    let accepted = advance_surface_with_accommodation(&fixture::phase(), &binding, solid,
        &gas, &catalog, |h| fixture::temperature(h, chart), &law, dt, fixture::controls(), || false).unwrap();
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
    // Independent ideal-diatomic flux-weighted energy is 3 R T, not bulk cv*T or cp*T.
    let expected = 0.5 * 100000.0 / (2.0 * core::f64::consts::PI * r * 300.0).sqrt()
        * (3.0 * r * 100.0) * solid.area_m2 * dt;
    fixture::near(accepted.reflection().heat_into_gas_j, expected, 1e-12);
    let end = accepted.exchange();
    assert_eq!(end.solid().mass_kg, solid.mass_kg);
    assert_eq!(end.gas().mass_kg(), gas.mass_kg());
    assert_eq!(end.step().transferred_mass_kg, 0.0);
    assert_eq!(end.step().recession_m, 0.0);
    fixture::near(end.step().surface_temperature_k[1], 400.0 - expected / (solid.mass_kg * chart.1), 1e-9);
    fixture::near(end.step().gas_temperature_k[1], 300.0 + expected / (gas.mass_kg() * 2.5 * r), 1e-8);
    fixture::near(end.step().solid_energy_change_j, -expected, 1e-8);
    fixture::near(end.step().gas_energy_change_j, expected, 1e-7);
    fixture::near(end.step().energy_residual_j, 0.0, 1e-7);
    assert_eq!(end.gas().point().axes().get("local-lot"), Some(&9.0));
}

#[test]
fn zero_thermal_accommodation_and_fully_sticking_vapor_retain_the_original_path() {
    let catalog = fixture::catalog(0.0);
    let chart = fixture::chart(&catalog, &fixture::binding(0.0));
    let gas = fixture::gas(&catalog, 350.0, 1000.0);
    let solid = fixture::solid(400.0, chart);
    for (sticking, vapor_heat) in [(0.1, 0.0), (1.0, 1.0)] {
        let mut binding = fixture::binding(0.0);
        binding.accommodation_coefficient = sticking;
        let law = accommodation(&catalog, 0.0, vapor_heat);
        let original = advance_surface_with_gas(&fixture::phase(), &binding, solid, &gas, &catalog,
            |h| fixture::temperature(h, chart), 0.001, fixture::controls(), || false).unwrap();
        let extended = advance_surface_with_accommodation(&fixture::phase(), &binding, solid, &gas,
            &catalog, |h| fixture::temperature(h, chart), &law, 0.001, fixture::controls(), || false).unwrap();
        assert_eq!(extended.reflection().heat_into_gas_j, 0.0);
        assert_eq!(extended.exchange(), &original);
    }
}

#[test]
fn reflected_vapor_uses_partial_pressure_and_only_the_nonsticking_fraction() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 300.0, 1000.0);
    let law = accommodation(&catalog, 0.0, 0.8);
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015;
    let mut results = Vec::new();
    for sticking in [0.25, 0.75] {
        let mut binding = fixture::binding(0.0);
        binding.accommodation_coefficient = sticking;
        let result = law.evaluate(&binding, &gas, 400.0, 1e-4, 0.01, || false).unwrap();
        let expected = 0.8 * (1.0 - sticking) * 1000.0
            / (2.0 * core::f64::consts::PI * r * 300.0).sqrt() * 3.5 * r * 100.0 * 1e-6;
        fixture::near(result.heat_into_gas_j, expected, 1e-12);
        assert_eq!(result.species[0].heat_into_gas_j, 0.0);
        assert_eq!(result.species[1].species, fixture::species());
        results.push(result.heat_into_gas_j);
    }
    fixture::near(results[0], 3.0 * results[1], 1e-12);
    let equilibrium = law.evaluate(&fixture::binding(0.0), &gas, 300.0, 1e-4, 0.01, || false).unwrap();
    assert_eq!(equilibrium.heat_into_gas_j, 0.0);
    assert!(law.evaluate(&fixture::binding(0.0), &gas, 250.0, 1e-4, 0.01, || false).unwrap().heat_into_gas_j < 0.0);
}

#[test]
fn variable_cp_and_large_formation_offsets_share_the_original_sensible_difference() {
    let id = SpeciesId::new("A").unwrap();
    for slope in [0.0, 0.001] {
        let mut results = Vec::new();
        for offset in [0.0, 1e12] {
            let catalog = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
                vec![species_with_energy_reference("A", 0.04, 2.5, slope, offset)],
                Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
            ).unwrap(), vec![transport("A", 1.7e-5, 0.025)], None).unwrap();
            let gas = ClosedGasParcel::new(catalog.clone(), point(300.0, 100000.0), 0.01).unwrap();
            let law = ThermalAccommodation::new(&catalog, &[(id.clone(), 0.6)], [200.0, 900.0]).unwrap();
            let t = 300.000001;
            let heat = law.evaluate(&fixture::binding(0.0), &gas, t, 1e-4, 0.001, || false).unwrap();
            let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.04;
            let de = r * (t - 300.0) * (2.0 + 0.5 * slope * (t + 300.0));
            let expected = 0.6 * 100000.0 / (2.0 * core::f64::consts::PI * r * 300.0).sqrt() * de * 1e-7;
            fixture::near(heat.heat_into_gas_j, expected, 1e-12 * expected.abs());
            results.push(heat.heat_into_gas_j);
        }
        assert_eq!(results[0], results[1]);
    }
}

#[test]
fn phase_and_reflected_heat_are_debited_once_at_one_endpoint() {
    let catalog = fixture::catalog(-1e7);
    let binding = fixture::binding(-1e7);
    let chart = fixture::chart(&catalog, &binding);
    let solid = fixture::solid(400.0, chart);
    let gas = fixture::gas(&catalog, 350.0, 5000.0);
    let law = accommodation(&catalog, 0.5, 0.6);
    let original = advance_surface_with_gas(&fixture::phase(), &binding, solid, &gas, &catalog,
        |h| fixture::temperature(h, chart), 0.001, fixture::controls(), || false).unwrap();
    let extended = advance_surface_with_accommodation(&fixture::phase(), &binding, solid, &gas,
        &catalog, |h| fixture::temperature(h, chart), &law, 0.001, fixture::controls(), || false).unwrap();
    let end = extended.exchange();
    assert_eq!(end.solid().mass_kg, original.solid().mass_kg);
    assert_eq!(end.step().transferred_mass_kg, original.step().transferred_mass_kg);
    fixture::near(end.step().energy_into_gas_j,
        original.step().energy_into_gas_j + extended.reflection().heat_into_gas_j, 1e-10);
    fixture::near(end.step().solid_energy_change_j + end.step().energy_into_gas_j, 0.0, 1e-7);
    fixture::near(end.step().energy_residual_j, 0.0, 1e-7);
    for id in [SpeciesId::new("A").unwrap(), fixture::species()] {
        fixture::near(end.gas().species_mass_kg(&id).unwrap(), original.gas().species_mass_kg(&id).unwrap(), 1e-14);
    }
}

#[test]
fn large_phase_datum_preserves_molecular_cooling_and_reflected_heat() {
    let base = fixture::catalog(0.0);
    let chart = fixture::chart(&base, &fixture::binding(0.0));
    let solid = fixture::solid(400.0, chart);
    let mut original: Option<(f64, f64, f64, f64, f64)> = None;
    for reference in [0.0, -1e14, 1e14] {
        let catalog = fixture::catalog(reference);
        let binding = fixture::binding(reference);
        let gas = fixture::gas(&catalog, 399.5, 5000.0);
        let saved = gas.clone();
        let law = accommodation(&catalog, 0.5, 0.6);
        let mut controls = fixture::controls();
        // Only the absolute source compatibility allowance reflects the datum's
        // ULP. The physical finite-energy and mass budgets remain unchanged.
        controls.latent_compatibility_tolerance_j_kg = 0.125;
        let result = advance_surface_with_accommodation(&fixture::phase(), &binding,
            solid, &gas, &catalog, |h| fixture::temperature(h, chart), &law,
            0.001, controls, || false).unwrap();
        let end = result.exchange();
        let step = end.step();
        let q = result.reflection().heat_into_gas_j;
        let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015;
        let ra = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
        let dm = step.transferred_mass_kg;
        let molecular_difference = 3.5 * r * (400.0 - 399.5);
        let thermal = dm * (r * 399.5 + 4.0 * r * 0.5 - 0.5 * r * 400.0)
            + step.incoming_mass_kg * molecular_difference + q;
        let final_capacity = gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap() * 2.5 * ra
            + (gas.species_mass_kg(&fixture::species()).unwrap() + dm) * 3.0 * r;
        fixture::near(step.gas_temperature_k[1], 399.5 + thermal / final_capacity, 2e-8);
        let alignment = (catalog.species_specific_enthalpy_j_kg(&fixture::species(), 400.0).unwrap()
            - (solid.specific_enthalpy_j_kg + binding.solid_enthalpy_offset_j_kg)) - 200000.0;
        let solid_heat = dm * (-200000.0 - alignment + 0.5 * r * 400.0)
            - step.incoming_mass_kg * molecular_difference - q;
        fixture::near(step.surface_temperature_k[1],
            400.0 + solid_heat / (end.solid().mass_kg * chart.1), 2e-9);
        assert!(q > 1e-5 && dm < 0.0);
        assert!(step.energy_residual_j.abs() <= controls.gas.energy_tolerance_j);
        let alignment_temperature = (dm * alignment).abs() / (end.solid().mass_kg * chart.1);
        let state = (step.gas_temperature_k[1], step.surface_temperature_k[1], dm, q, alignment_temperature);
        if let Some(first) = original {
            fixture::near(state.0, first.0, 2e-8);
            fixture::near(state.1, first.1, 2e-9 + state.4 + first.4);
            assert_eq!((state.2, state.3), (first.2, first.3));
        } else { original = Some(state); }
        assert_eq!(gas, saved);
    }
}

#[test]
fn reflected_finite_store_trajectory_refines_toward_independent_analytic_solution() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let chart = fixture::chart(&catalog, &binding);
    let initial_solid = fixture::solid(400.0, chart);
    let initial_gas = fixture::gas(&catalog, 300.0, 0.0);
    let law = accommodation(&catalog, 0.5, 0.0);
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
    let cs = initial_solid.mass_kg * chart.1;
    let cg = initial_gas.mass_kg() * 2.5 * r;
    let equilibrium = (cs * 400.0 + cg * 300.0) / (cs + cg);
    let k = 0.5 * initial_solid.area_m2 * initial_gas.mass_kg() / 0.01
        * (r / (2.0 * core::f64::consts::PI)).sqrt() * 3.0 * r;
    let x = (300.0 / equilibrium).sqrt();
    let z = 0.5 * fs_math::det::ln((1.0 + x) / (1.0 - x))
        + k * (1.0 + cg / cs) * equilibrium.sqrt() * 0.1 / (2.0 * cg);
    let v = fs_math::det::tanh(z);
    let exact_temperature = equilibrium * v * v;
    let mut errors = Vec::new();
    for count in [4, 8, 16] {
        let mut solid = initial_solid;
        let mut gas = initial_gas.clone();
        for _ in 0..count {
            let candidate = advance_surface_with_accommodation(&fixture::phase(), &binding,
                solid, &gas, &catalog, |h| fixture::temperature(h, chart), &law,
                0.1 / f64::from(count), fixture::controls(), || false).unwrap();
            let (exchange, _) = candidate.into_parts();
            let (next_solid, next_gas, _, _) = exchange.into_parts();
            solid = next_solid;
            gas = next_gas;
        }
        errors.push((gas.state().state().temperature - exact_temperature).abs());
    }
    assert!(errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1], "{errors:?}");
}

#[test]
fn missing_coefficients_invalid_domains_and_late_cancellation_preserve_inputs() {
    let catalog = fixture::catalog(0.0);
    let a = SpeciesId::new("A").unwrap();
    let b = fixture::species();
    assert!(ThermalAccommodation::new(&catalog, &[(a.clone(), 0.5)], [200.0, 900.0]).is_err());
    assert!(ThermalAccommodation::new(&catalog, &[(a.clone(), 0.5), (a.clone(), 0.5)], [200.0, 900.0]).is_err());
    assert!(ThermalAccommodation::new(&catalog, &[(a.clone(), 0.5), (b.clone(), 1.1)], [200.0, 900.0]).is_err());
    assert!(ThermalAccommodation::new(&catalog, &[(a.clone(), 0.5), (b.clone(), 0.5)], [900.0, 200.0]).is_err());
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let solid = fixture::solid(400.0, chart);
    let gas = fixture::gas(&catalog, 350.0, 1000.0);
    let saved = gas.clone();
    let law = accommodation(&catalog, 0.5, 0.5);
    let restricted = ThermalAccommodation::new(&catalog, &[(a.clone(), 0.5), (b.clone(), 0.5)], [350.0, 400.0]).unwrap();
    assert!(restricted.evaluate(&binding, &gas, 401.0, 1e-4, 0.01, || false).is_err());
    let reversed = ThermalAccommodation::new(&catalog, &[(a, 0.5), (b, 0.5)], [200.0, 900.0]).unwrap();
    assert_eq!(law, reversed);
    let mut calls = 0;
    advance_surface_with_accommodation(&fixture::phase(), &binding, solid, &gas, &catalog,
        |h| fixture::temperature(h, chart), &law, 0.001, fixture::controls(), || { calls += 1; false }).unwrap();
    for stop in [1, calls / 2, calls] {
        let mut current = 0;
        assert!(advance_surface_with_accommodation(&fixture::phase(), &binding, solid, &gas,
            &catalog, |h| fixture::temperature(h, chart), &law, 0.001, fixture::controls(),
            || { current += 1; current == stop }).is_err());
        assert_eq!(gas, saved);
        assert_eq!(solid, fixture::solid(400.0, chart));
    }
}
