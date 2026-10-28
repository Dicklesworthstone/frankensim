//! Joint exchange through the original gas/material owners and compatible charts.
use super::*;
use crate::sublimation::gas::{advance_surface_with_gas, tests as fixture};
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

fn controls() -> GasChamberSolve {
    GasChamberSolve {
        surface: fixture::controls(),
        maximum_surfaces: 32,
    }
}
fn surface<'a>(
    id: u64,
    law: &'a dyn SublimationLaw,
    binding: &'a GasSurfaceBinding,
    solid: SolidSurfaceState,
    caloric: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
) -> GasChamberSurface<'a> {
    GasChamberSurface {
        id,
        law,
        binding,
        solid,
        temperature_from_enthalpy: caloric,
    }
}

#[test]
fn hot_and_cold_surfaces_share_one_evolving_gas_and_are_permutation_invariant() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let warm = fixture::solid(500.0, chart);
    let cold = fixture::solid(300.0, chart);
    let rows = [
        surface(20, &law, &binding, warm, &caloric),
        surface(10, &law, &binding, cold, &caloric),
    ];
    let accepted =
        advance_chamber_with_gas(&rows, &gas, &catalog, 0.1, controls(), || false).unwrap();
    let permuted =
        advance_chamber_with_gas(&[rows[1], rows[0]], &gas, &catalog, 0.1, controls(), || {
            false
        })
        .unwrap();
    assert_eq!(accepted, permuted);
    assert_eq!(
        accepted
            .surfaces()
            .iter()
            .map(|row| row.id)
            .collect::<Vec<_>>(),
        [10, 20]
    );
    assert!(accepted.surfaces()[0].transferred_mass_kg < 0.0);
    assert!(accepted.surfaces()[1].transferred_mass_kg > 0.0);
    assert!(accepted.surfaces()[0].solid.mass_kg > cold.mass_kg);
    assert!(accepted.surfaces()[1].solid.mass_kg < warm.mass_kg);
    assert_ne!(
        accepted.gas().state().state().temperature,
        gas.state().state().temperature
    );
    let carrier = SpeciesId::new("A").unwrap();
    fixture::near(
        accepted.gas().species_mass_kg(&carrier).unwrap(),
        gas.species_mass_kg(&carrier).unwrap(),
        1e-14,
    );
    fixture::near(accepted.step().energy_residual_j, 0.0, 1e-7);
    fixture::near(accepted.step().mass_residual_kg, 0.0, 1e-14);
    assert_eq!(accepted.step().species.len(), 1);
    fixture::near(accepted.step().species[0].mass_residual_kg, 0.0, 1e-14);
    for surface in accepted.surfaces() {
        fixture::near(surface.energy_residual_j, 0.0, 1e-7);
    }
}

#[test]
fn large_datum_counterstreams_and_carrier_reflection_share_one_physical_target() {
    let base = fixture::catalog(0.0);
    let chart = fixture::chart(&base, &fixture::binding(0.0));
    let caloric = |h| fixture::temperature(h, chart);
    let phase = fixture::phase();
    let r = super::super::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015;
    let ra = super::super::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
    let mut original: Option<(f64, f64, Vec<(f64, f64)>)> = None;
    for reference in [0.0, -1e14, 1e14] {
        let catalog = fixture::catalog(reference);
        let binding = fixture::binding(reference);
        let gas = fixture::gas(&catalog, 400.0, 1000.0);
        let rows = [surface(2, &phase, &binding, fixture::solid(400.5, chart), &caloric),
            surface(1, &phase, &binding, fixture::solid(399.75, chart), &caloric)];
        let accommodation = ThermalAccommodation::new(&catalog, &[
            (SpeciesId::new("A").unwrap(), 0.5), (fixture::species(), 0.6),
        ], [200.0, 900.0]).unwrap();
        let maps = BTreeMap::from([(1, accommodation.clone()), (2, accommodation)]);
        let mut controls = controls();
        controls.surface.latent_compatibility_tolerance_j_kg = 0.125;
        let (result, reflection) = advance_chamber_with_accommodation(&rows, &maps,
            &gas, &catalog, 0.001, controls, || false).unwrap();
        let reversed = advance_chamber_with_accommodation(&[rows[1], rows[0]], &maps,
            &gas, &catalog, 0.001, controls, || false).unwrap();
        assert_eq!((&result, &reflection), (&reversed.0, &reversed.1));
        assert!(result.surfaces()[0].transferred_mass_kg < 0.0);
        assert!(result.surfaces()[1].transferred_mass_kg > 0.0);
        let mut thermal = 0.0;
        let mut solid_temperatures = Vec::new();
        for row in result.surfaces() {
            let ts = row.surface_temperature_k[0];
            thermal += row.transferred_mass_kg * (r * 400.0 + 4.0 * r * (ts - 400.0) - 0.5 * r * ts)
                + row.incoming_mass_kg * 3.5 * r * (ts - 400.0)
                + reflection[&row.id].heat_into_gas_j;
            let initial = rows.iter().find(|source| source.id == row.id).unwrap().solid;
            let alignment = (catalog.species_specific_enthalpy_j_kg(&fixture::species(), ts).unwrap()
                - (initial.specific_enthalpy_j_kg + binding.solid_enthalpy_offset_j_kg)) - 200000.0;
            let solid_heat = row.transferred_mass_kg * (-200000.0 - alignment + 0.5 * r * ts)
                - row.incoming_mass_kg * 3.5 * r * (ts - 400.0) - reflection[&row.id].heat_into_gas_j;
            fixture::near(row.surface_temperature_k[1],
                ts + solid_heat / (row.solid.mass_kg * chart.1), 2e-9);
            solid_temperatures.push((row.surface_temperature_k[1],
                (row.transferred_mass_kg * alignment).abs() / (row.solid.mass_kg * chart.1)));
            assert!(row.energy_residual_j.abs() <= controls.surface.gas.energy_tolerance_j);
        }
        let capacity = gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap() * 2.5 * ra
            + (gas.species_mass_kg(&fixture::species()).unwrap() + result.step().transferred_mass_kg) * 3.0 * r;
        fixture::near(result.gas().state().state().temperature, 400.0 + thermal / capacity, 2e-8);
        assert!(result.step().energy_residual_j.abs() <= controls.surface.gas.energy_tolerance_j);
        assert!(result.step().mass_residual_kg.abs() <= controls.surface.gas.mass_tolerance_kg);
        let state = (result.gas().state().state().temperature, result.step().transferred_mass_kg,
            solid_temperatures);
        if let Some(first) = &original {
            fixture::near(state.0, first.0, 2e-8);
            assert_eq!(state.1, first.1);
            for (a, b) in state.2.iter().zip(&first.2) { fixture::near(a.0, b.0, 2e-9 + a.1 + b.1); }
        } else { original = Some(state); }
    }
}

#[test]
fn heterogeneous_species_and_caloric_sources_preserve_each_carrier_and_species() {
    let b = SpeciesId::new("B").unwrap();
    let c = SpeciesId::new("C").unwrap();
    let catalog = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.02897, 3.5, 0.0, 0.0),
                species_with_energy_reference("B", 0.018015, 4.0, 0.0, 5e6),
                species_with_energy_reference("C", 0.044, 3.5, 0.0, -1e7),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.98, 0.01, 0.01]).unwrap(),
        )
        .unwrap(),
        vec![
            transport("A", 1.7e-5, 0.025),
            transport("B", 2e-5, 0.03),
            transport("C", 1.5e-5, 0.02),
        ],
        Some(&b),
    )
    .unwrap();
    let gas = ClosedGasParcel::new(catalog.clone(), point(400.0, 100000.0), 0.01).unwrap();
    let bindings = [
        GasSurfaceBinding {
            species: b.clone(),
            accommodation_coefficient: 0.1,
            solid_enthalpy_offset_j_kg: 5e6,
        },
        GasSurfaceBinding {
            species: c.clone(),
            accommodation_coefficient: 0.1,
            solid_enthalpy_offset_j_kg: -1e7,
        },
    ];
    let laws = [
        fixture::phase(),
        fixture::Phase {
            molar_mass: 0.044,
            reference_temperature_k: 400.0,
        },
    ];
    let charts = [
        fixture::chart(&catalog, &bindings[0]),
        fixture::chart(&catalog, &bindings[1]),
    ];
    let caloric_b = |h| fixture::temperature(h, charts[0]);
    let caloric_c = |h| fixture::temperature(h, charts[1]);
    let rows = [
        surface(
            2,
            &laws[1],
            &bindings[1],
            fixture::solid(300.0, charts[1]),
            &caloric_c,
        ),
        surface(
            1,
            &laws[0],
            &bindings[0],
            fixture::solid(500.0, charts[0]),
            &caloric_b,
        ),
    ];
    let result =
        advance_chamber_with_gas(&rows, &gas, &catalog, 0.1, controls(), || false).unwrap();
    assert!(result.surfaces()[0].transferred_mass_kg > 0.0);
    assert!(result.surfaces()[1].transferred_mass_kg < 0.0);
    assert_eq!(
        result
            .step()
            .species
            .iter()
            .map(|row| row.species.clone())
            .collect::<Vec<_>>(),
        [b, c]
    );
    let carrier = SpeciesId::new("A").unwrap();
    fixture::near(
        result.gas().species_mass_kg(&carrier).unwrap(),
        gas.species_mass_kg(&carrier).unwrap(),
        1e-14,
    );
    for row in &result.step().species {
        fixture::near(
            row.vapor_mass_kg[1] - row.vapor_mass_kg[0],
            row.transferred_mass_kg,
            1e-14,
        );
        fixture::near(row.mass_residual_kg, 0.0, 1e-14);
    }
    assert_eq!(
        result.gas().model().thermodynamics().components(),
        catalog.thermodynamics().components()
    );
    fixture::near(result.step().energy_residual_j, 0.0, 1e-7);
}

#[test]
fn aggregate_gross_withdrawal_cannot_spend_the_same_original_vapor_twice() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let saved = gas.clone();
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let solid = fixture::solid(300.0, chart);
    // Each surface separately fits the 25% allowance. Their simultaneous gross
    // withdrawal exceeds it even though net deposition consumes much less.
    advance_surface_with_gas(
        &law,
        &binding,
        solid,
        &gas,
        &catalog,
        caloric,
        1.0,
        controls().surface,
        || false,
    )
    .unwrap();
    let rows = [
        surface(1, &law, &binding, solid, &caloric),
        surface(2, &law, &binding, solid, &caloric),
    ];
    let refused = advance_chamber_with_gas(&rows, &gas, &catalog, 1.0, controls(), || false);
    let Err(GasSurfaceError::InventoryLimit { maximum_duration_s }) = refused else {
        panic!("combined original-vapor budget must refuse");
    };
    assert!(maximum_duration_s > 0.0 && maximum_duration_s < 1.0);
    assert_eq!(gas, saved);
    let accepted = advance_chamber_with_gas(
        &rows,
        &gas,
        &catalog,
        0.99 * maximum_duration_s,
        controls(),
        || false,
    )
    .unwrap();
    assert!(
        accepted.step().species[0].incoming_mass_kg / accepted.step().species[0].vapor_mass_kg[0]
            <= controls().surface.maximum_gross_inventory_fraction
    );
}

#[test]
fn one_surface_is_the_original_single_surface_update_exactly() {
    let catalog = fixture::catalog(-1e7);
    let binding = fixture::binding(-1e7);
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let law = fixture::phase();
    let solid = fixture::solid(400.0, chart);
    for pressure in [0.0, 500.0, 5000.0] {
        let gas = fixture::gas(&catalog, 400.0, pressure);
        let one = advance_surface_with_gas(
            &law,
            &binding,
            solid,
            &gas,
            &catalog,
            caloric,
            0.1,
            controls().surface,
            || false,
        )
        .unwrap();
        let joint = advance_chamber_with_gas(
            &[surface(7, &law, &binding, solid, &caloric)],
            &gas,
            &catalog,
            0.1,
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(one.gas_exchange(), joint.gas_exchange());
        assert_eq!(one.solid(), joint.surfaces()[0].solid);
        assert_eq!(
            one.step().transferred_mass_kg,
            joint.surfaces()[0].transferred_mass_kg
        );
        assert_eq!(one.step().energy_into_gas_j, joint.step().energy_into_gas_j);
        assert_eq!(one.step().energy_residual_j, joint.step().energy_residual_j);
        assert_eq!(one.step().mass_residual_kg, joint.step().mass_residual_kg);
    }
}

#[test]
fn common_formation_reference_survives_opposing_surface_energy_aggregation() {
    let base = fixture::catalog(0.0);
    let shifted = fixture::catalog(1e8);
    let law = fixture::phase();
    let run = |catalog: &ThermochemicalGasModel, offset: f64| {
        let binding = fixture::binding(offset);
        let chart = fixture::chart(catalog, &binding);
        let caloric = |h| fixture::temperature(h, chart);
        let gas = fixture::gas(catalog, 400.0, 1000.0);
        let rows = [
            surface(1, &law, &binding, fixture::solid(500.0, chart), &caloric),
            surface(2, &law, &binding, fixture::solid(300.0, chart), &caloric),
        ];
        let result =
            advance_chamber_with_gas(&rows, &gas, catalog, 0.1, controls(), || false).unwrap();
        (
            result.surfaces().to_vec(),
            result.gas().state().state().temperature,
            result.step().clone(),
        )
    };
    let a = run(&base, 0.0);
    let b = run(&shifted, 1e8);
    fixture::near(a.1, b.1, 1e-7);
    for (left, right) in a.0.iter().zip(&b.0) {
        fixture::near(
            left.surface_temperature_k[1],
            right.surface_temperature_k[1],
            1e-7,
        );
        assert_eq!(left.transferred_mass_kg, right.transferred_mass_kg);
    }
    fixture::near(
        b.2.energy_into_gas_j - a.2.energy_into_gas_j,
        1e8 * a.2.transferred_mass_kg,
        1e-7,
    );
    fixture::near(a.2.energy_residual_j, 0.0, 1e-7);
    fixture::near(b.2.energy_residual_j, 0.0, 1e-7);
}

#[test]
fn duplicate_ids_budgets_late_caloric_failure_and_cancellation_publish_nothing() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let saved = gas.clone();
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let caloric = |h| fixture::temperature(h, chart);
    let warm = fixture::solid(500.0, chart);
    let cold = fixture::solid(300.0, chart);
    let rows = [
        surface(1, &law, &binding, warm, &caloric),
        surface(2, &law, &binding, cold, &caloric),
    ];
    assert!(
        advance_chamber_with_gas(&[rows[0], rows[0]], &gas, &catalog, 0.1, controls(), || {
            false
        })
        .is_err()
    );
    assert!(
        advance_chamber_with_gas(
            &rows,
            &gas,
            &catalog,
            0.1,
            GasChamberSolve {
                maximum_surfaces: 1,
                ..controls()
            },
            || false
        )
        .is_err()
    );
    let broken = |h| {
        fixture::temperature(h, chart).map(|t| {
            if h == cold.specific_enthalpy_j_kg {
                t
            } else {
                t + 0.5
            }
        })
    };
    let bad = [rows[0], surface(2, &law, &binding, cold, &broken)];
    assert!(matches!(
        advance_chamber_with_gas(&bad, &gas, &catalog, 0.1, controls(), || false),
        Err(GasSurfaceError::CaloricMismatch { .. })
    ));
    let empty = advance_chamber_with_gas(&[], &gas, &catalog, 0.1, controls(), || false).unwrap();
    assert_eq!(empty.gas(), &gas);
    assert!(empty.surfaces().is_empty());
    assert_eq!(empty.step().energy_residual_j, 0.0);
    let mut calls = 0;
    advance_chamber_with_gas(&rows, &gas, &catalog, 0.1, controls(), || {
        calls += 1;
        false
    })
    .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(
            advance_chamber_with_gas(&rows, &gas, &catalog, 0.1, controls(), || {
                count += 1;
                count == stop
            })
            .is_err()
        );
        assert_eq!(gas, saved);
        assert_eq!(rows[0].solid, warm);
        assert_eq!(rows[1].solid, cold);
    }
}
