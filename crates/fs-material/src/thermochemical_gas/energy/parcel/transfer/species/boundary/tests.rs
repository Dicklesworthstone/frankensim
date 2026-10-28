//! Synthetic caloric fixtures exercise the original thermochemical/material owners.
use super::*;
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis};

fn controls() -> GasSpeciesBoundarySolve {
    GasSpeciesBoundarySolve {
        gas: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1000.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        mass_tolerance_kg: 1e-14,
        energy_tolerance_j: 1e-7,
    }
}
fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn catalog(reference: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.02897, 3.5, 0.0, 0.0),
                species_with_energy_reference("B", 0.018015, 4.0, 0.001, reference),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.9, 0.1]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.03)],
        Some(&id("B")),
    )
    .unwrap()
}
fn parcel(model: &ThermochemicalGasModel) -> ClosedGasParcel {
    ClosedGasParcel::new(
        model.clone(),
        point(400.0, 100000.0).with("local-lot", 7.0).unwrap(),
        0.001,
    )
    .unwrap()
}
fn stream(
    model: &ThermochemicalGasModel,
    name: &str,
    mass: f64,
    temperature: f64,
) -> SpeciesBoundaryTransfer {
    let species = id(name);
    let enthalpy = model
        .species_specific_enthalpy_j_kg(&species, temperature)
        .unwrap();
    SpeciesBoundaryTransfer {
        species,
        mass_into_gas_kg: mass,
        energy_into_gas_j: mass * enthalpy,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn influx_and_deposition_use_original_partial_enthalpy_and_preserve_carrier() {
    let catalog = catalog(0.0);
    let gas = parcel(&catalog);
    let vapor = gas.species_mass_kg(&id("B")).unwrap();
    let carrier = gas.species_mass_kg(&id("A")).unwrap();
    near(
        gas.species_partial_pressure_pa(&id("B")).unwrap(),
        10000.0,
        1e-9,
    );
    for (fraction, temperature) in [(0.2, 600.0), (-0.2, 400.0)] {
        let row = stream(&catalog, "B", fraction * vapor, temperature);
        let candidate = gas
            .exchange_species_boundary(&catalog, &[row.clone()], controls(), || false)
            .unwrap();
        let next = candidate.gas();
        near(next.species_mass_kg(&id("A")).unwrap(), carrier, 1e-14);
        near(
            next.species_mass_kg(&id("B")).unwrap(),
            vapor + row.mass_into_gas_kg,
            1e-14,
        );
        near(
            candidate.step().stored_energy_change_j,
            row.energy_into_gas_j,
            1e-7,
        );
        assert!(candidate.step().maximum_species_mass_residual_kg <= 1e-14);
        assert_eq!(next.volume_m3(), gas.volume_m3());
        assert_eq!(next.point().axes().get("local-lot"), Some(&7.0));
        assert_eq!(next.model().water_species(), Some(&id("B")));
        assert_eq!(
            next.model().transport_materials(),
            catalog.transport_materials()
        );
        assert_ne!(
            next.state().state().temperature,
            gas.state().state().temperature
        );
        assert_ne!(next.state().state().pressure, gas.state().state().pressure);
        if fraction > 0.0 {
            assert!(next.state().state().temperature > gas.state().state().temperature);
        } else {
            assert!(next.state().state().temperature < gas.state().state().temperature);
        }
    }
}

#[test]
fn exact_vapor_depletion_and_reentry_keep_the_complete_original_catalog() {
    let catalog = catalog(-1.0e7);
    let gas = parcel(&catalog);
    let vapor = gas.species_mass_kg(&id("B")).unwrap();
    let removal = stream(&catalog, "B", -vapor, 400.0);
    let dry = gas
        .exchange_species_boundary(&catalog, &[removal.clone()], controls(), || false)
        .unwrap();
    assert!(std::ptr::eq(dry.catalog(), &catalog));
    assert_eq!(dry.catalog().thermodynamics().components().len(), 2);
    assert_eq!(dry.gas().model().thermodynamics().components().len(), 1);
    assert_eq!(dry.gas().model().water_species(), None);
    assert_eq!(dry.gas().species_mass_kg(&id("B")).unwrap(), 0.0);
    assert_eq!(
        dry.gas().species_partial_pressure_pa(&id("B")).unwrap(),
        0.0
    );
    let replacement = SpeciesBoundaryTransfer {
        species: id("B"),
        mass_into_gas_kg: vapor,
        energy_into_gas_j: -removal.energy_into_gas_j,
    };
    let wet = dry
        .gas()
        .exchange_species_boundary(dry.catalog(), &[replacement], controls(), || false)
        .unwrap();
    assert_eq!(
        wet.gas().model().thermodynamics().components(),
        catalog.thermodynamics().components()
    );
    assert_eq!(wet.gas().model().water_species(), Some(&id("B")));
    near(wet.gas().mass_kg(), gas.mass_kg(), 1e-14);
    near(wet.gas().state().state().temperature, 400.0, 1e-7);
    near(wet.gas().state().state().pressure, 100000.0, 1e-4);
    near(
        dry.step().stored_energy_change_j + wet.step().stored_energy_change_j,
        0.0,
        2e-7,
    );
}

#[test]
fn species_energy_reference_changes_energy_accounts_without_changing_physics() {
    let base = catalog(0.0);
    let shifted = catalog(-1.0e8);
    let original = [parcel(&base), parcel(&shifted)];
    let mass = 0.1 * original[0].species_mass_kg(&id("B")).unwrap();
    let trials = [&base, &shifted]
        .into_iter()
        .enumerate()
        .map(|(index, source)| {
            original[index]
                .exchange_species_boundary(
                    source,
                    &[stream(source, "B", mass, 550.0)],
                    controls(),
                    || false,
                )
                .unwrap()
        })
        .collect::<Vec<_>>();
    near(
        trials[0].gas().state().state().temperature,
        trials[1].gas().state().state().temperature,
        1e-7,
    );
    near(
        trials[0].gas().state().state().pressure,
        trials[1].gas().state().state().pressure,
        1e-4,
    );
    near(
        trials[1].step().energy_into_gas_j - trials[0].step().energy_into_gas_j,
        -1.0e8 * mass,
        1e-8,
    );
    near(
        shifted
            .species_specific_enthalpy_j_kg(&id("B"), 550.0)
            .unwrap()
            - base
                .species_specific_enthalpy_j_kg(&id("B"), 550.0)
                .unwrap(),
        -1.0e8,
        1e-7,
    );
}

#[test]
fn zero_net_mass_counterflow_matches_original_heat_inversion() {
    let catalog = catalog(-1.0e7);
    let gas = parcel(&catalog);
    let row = SpeciesBoundaryTransfer {
        species: id("B"),
        mass_into_gas_kg: 0.0,
        energy_into_gas_j: 3.0,
    };
    let trial = gas
        .exchange_species_boundary(&catalog, &[row], controls(), || false)
        .unwrap();
    let mut reference = gas.clone();
    reference.add_heat(3.0, controls().gas, || false).unwrap();
    near(
        trial.gas().state().state().temperature,
        reference.state().state().temperature,
        1e-8,
    );
    near(
        trial.gas().state().state().pressure,
        reference.state().state().pressure,
        1e-5,
    );
    assert_eq!(trial.gas().mass_kg(), gas.mass_kg());
    near(trial.step().stored_energy_change_j, 3.0, 1e-7);
}

#[test]
fn canonical_rows_and_zero_transfer_preserve_the_original_owner() {
    let catalog = catalog(0.0);
    let gas = parcel(&catalog);
    let rows = [
        stream(&catalog, "A", -1e-6, 400.0),
        stream(&catalog, "B", 2e-6, 500.0),
    ];
    let first = gas
        .exchange_species_boundary(&catalog, &rows, controls(), || false)
        .unwrap();
    let second = gas
        .exchange_species_boundary(
            &catalog,
            &[rows[1].clone(), rows[0].clone()],
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(first, second);
    let empty = gas
        .exchange_species_boundary(&catalog, &[], controls(), || false)
        .unwrap();
    assert_eq!(empty.gas(), &gas);
    assert_eq!(empty.step().iterations, 0);
}

#[test]
fn invalid_catalog_inventory_energy_and_cancelled_trials_publish_nothing() {
    let catalog = catalog(0.0);
    let gas = parcel(&catalog);
    let saved = gas.clone();
    let valid = stream(&catalog, "B", 1e-6, 500.0);
    let vapor = gas.species_mass_kg(&id("B")).unwrap();
    let bad = [
        vec![SpeciesBoundaryTransfer {
            energy_into_gas_j: f64::NAN,
            ..valid.clone()
        }],
        vec![SpeciesBoundaryTransfer {
            mass_into_gas_kg: f64::INFINITY,
            ..valid.clone()
        }],
        vec![SpeciesBoundaryTransfer {
            species: id("foreign"),
            ..valid.clone()
        }],
        vec![valid.clone(), valid.clone()],
        vec![stream(&catalog, "B", -2.0 * vapor, 400.0)],
        vec![SpeciesBoundaryTransfer {
            energy_into_gas_j: 1.0e9,
            ..valid.clone()
        }],
    ];
    for rows in bad {
        assert!(
            gas.exchange_species_boundary(&catalog, &rows, controls(), || false)
                .is_err()
        );
        assert_eq!(gas, saved);
    }
    let changed = self::catalog(1.0);
    assert!(
        gas.exchange_species_boundary(&changed, &[valid.clone()], controls(), || false)
            .is_err()
    );
    let dry_catalog = catalog
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap())
        .unwrap();
    assert!(
        gas.exchange_species_boundary(&dry_catalog, &[valid.clone()], controls(), || false)
            .is_err()
    );
    for temperature in [199.0, 2001.0, 0.0, f64::NAN] {
        assert!(
            catalog
                .species_specific_enthalpy_j_kg(&id("B"), temperature)
                .is_err()
        );
    }
    assert!(
        catalog
            .species_specific_enthalpy_j_kg(&id("foreign"), 400.0)
            .is_err()
    );
    let mut calls = 0;
    gas.exchange_species_boundary(&catalog, &[valid.clone()], controls(), || {
        calls += 1;
        false
    })
    .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert_eq!(
            gas.exchange_species_boundary(&catalog, &[valid.clone()], controls(), || {
                count += 1;
                count == stop
            }),
            Err(GasEnergyError::Cancelled)
        );
        assert_eq!(gas, saved);
    }
}

fn relative_stream(
    model: &ThermochemicalGasModel,
    name: &str,
    mass: f64,
    original_temperature: f64,
    upstream_temperature: f64,
) -> SpeciesBoundaryRelativeTransfer {
    let species = id(name);
    let change = model.species_specific_enthalpy_change_j_kg(
        &species, original_temperature, upstream_temperature,
    ).unwrap();
    SpeciesBoundaryRelativeTransfer {
        species, mass_into_gas_kg: mass, relative_energy_into_gas_j: mass * change,
    }
}

#[test]
fn relative_stream_matches_analytic_rigid_volume_heating_and_blowdown() {
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
    for reference in [0.0, -1.0e18, 1.0e18] {
        let model = ThermochemicalGasModel::new(
            FrozenIdealGasMixtureModelV1::new(
                vec![species_with_energy_reference("A", 0.02897, 3.5, 0.0, reference)],
                Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
            ).unwrap(), vec![transport("A", 1.7e-5, 0.025)], None,
        ).unwrap();
        let gas = parcel(&model);
        for (fraction, upstream, heat) in [(0.2, 650.0, 1.0e-4), (-0.2, 400.0, 0.0)] {
            let dm = fraction * gas.mass_kg();
            let mut row = relative_stream(&model, "A", dm, 400.0, upstream);
            row.relative_energy_into_gas_j += heat;
            // First law at fixed volume: final m*cv*dT = dm*R*T0 + dm*cp*dTin + Q.
            let expected = 400.0 + (dm * r * 400.0 + dm * 3.5 * r * (upstream - 400.0)
                + heat) / ((gas.mass_kg() + dm) * 2.5 * r);
            let next = gas.exchange_species_boundary_relative(
                &model, &[row], controls(), || false,
            ).unwrap();
            near(next.gas().state().state().temperature, expected, 2e-8);
            near(next.gas().state().state().pressure,
                (gas.mass_kg() + dm) * r * expected / gas.volume_m3(), 2e-5);
            near(next.step().energy_residual_j, 0.0, 1e-7);
            assert_eq!(next.gas().volume_m3(), gas.volume_m3());
            assert_eq!(next.gas().point().axes().get("local-lot"), Some(&7.0));
            assert_eq!(expected > 400.0, fraction > 0.0);
            if reference != 0.0 && heat != 0.0 {
                let absolute = dm * model.species_specific_enthalpy_j_kg(&id("A"), upstream).unwrap();
                assert_eq!(absolute + heat, absolute, "the thermal correction is below the absolute datum's ULP");
            }
        }
    }
}

#[test]
fn relative_mixed_streams_integrate_the_final_nonlinear_caloric_law() {
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
    let mut previous = None;
    for reference in [0.0, -1.0e18, 1.0e18] {
        let model = catalog(reference);
        let gas = parcel(&model);
        let ma = gas.species_mass_kg(&id("A")).unwrap();
        let mb = gas.species_mass_kg(&id("B")).unwrap();
        let da = -0.05 * ma;
        let db = 0.3 * mb;
        let mut rows = [relative_stream(&model, "A", da, 400.0, 400.0),
            relative_stream(&model, "B", db, 400.0, 650.0)];
        rows[1].relative_energy_into_gas_j += 0.002;
        let sensible = r * 400.0 * (da / 0.02897 + db / 0.018015)
            + db * r / 0.018015 * (4.0 * 250.0 + 0.0005 * (650.0_f64.powi(2) - 400.0_f64.powi(2)))
            + 0.002;
        // Integrate cv(T) analytically in the final mixture, then solve its quadratic.
        let constant = (ma + da) / 0.02897 * 2.5 + (mb + db) / 0.018015 * 3.0;
        let slope = (mb + db) / 0.018015 * 0.001;
        let at_original = constant + 400.0 * slope;
        let expected = 400.0 + 2.0 * sensible / r
            / (at_original + (at_original * at_original + 2.0 * slope * sensible / r).sqrt());
        let next = gas.exchange_species_boundary_relative(&model, &rows, controls(), || false).unwrap();
        let reversed = gas.exchange_species_boundary_relative(
            &model, &[rows[1].clone(), rows[0].clone()], controls(), || false,
        ).unwrap();
        assert_eq!(next, reversed);
        near(next.gas().state().state().temperature, expected, 2e-8);
        near(next.gas().species_mass_kg(&id("A")).unwrap(), ma + da, 1e-14);
        near(next.gas().species_mass_kg(&id("B")).unwrap(), mb + db, 1e-14);
        near(next.step().energy_residual_j, 0.0, 1e-7);
        if let Some((temperature, pressure)) = previous {
            near(next.gas().state().state().temperature, temperature, 1e-8);
            near(next.gas().state().state().pressure, pressure, 1e-5);
        }
        previous = Some((next.gas().state().state().temperature, next.gas().state().state().pressure));
    }
}

#[test]
fn relative_zero_mass_counterflow_preserves_small_heat_and_exact_noop() {
    let model = catalog(-1.0e18);
    let gas = parcel(&model);
    for heat in [-0.001, 0.0, 0.001] {
        let row = SpeciesBoundaryRelativeTransfer { species: id("B"),
            mass_into_gas_kg: 0.0, relative_energy_into_gas_j: heat };
        let next = gas.exchange_species_boundary_relative(&model, &[row], controls(), || false).unwrap();
        let mut heated = gas.clone();
        heated.add_heat(heat, controls().gas, || false).unwrap();
        near(next.gas().state().state().temperature, heated.state().state().temperature, 1e-8);
        assert_eq!(next.gas().mass_kg(), gas.mass_kg());
        assert_eq!(next.gas().model(), gas.model());
        near(next.step().stored_energy_change_j, heat, 1e-7);
        near(next.step().energy_residual_j, 0.0, 1e-7);
        if heat == 0.0 {
            assert_eq!(next.gas(), &gas);
            assert!(next.step().species.is_empty());
            assert_eq!(next.step().iterations, 0);
        }
    }
}

#[test]
fn relative_exact_depletion_and_reentry_keep_sources_at_large_reference() {
    for reference in [-1.0e18, 1.0e18] {
        let model = catalog(reference);
        let gas = parcel(&model);
        let vapor = gas.species_mass_kg(&id("B")).unwrap();
        let carrier = gas.species_mass_kg(&id("A")).unwrap();
        let dry = gas.exchange_species_boundary_relative(&model,
            &[relative_stream(&model, "B", -vapor, 400.0, 400.0)], controls(), || false).unwrap();
        assert!(std::ptr::eq(dry.catalog(), &model));
        assert_eq!(dry.catalog().thermodynamics().components().len(), 2);
        assert_eq!(dry.gas().model().thermodynamics().components().len(), 1);
        assert_eq!(dry.gas().species_mass_kg(&id("B")).unwrap(), 0.0);
        assert_eq!(dry.gas().species_partial_pressure_pa(&id("B")).unwrap(), 0.0);
        assert_eq!(dry.gas().model().water_species(), None);
        assert!(dry.gas().state().state().temperature < 400.0);
        // Restore precisely the removed upstream stream, using the NEW gas as reference.
        let replacement = relative_stream(dry.catalog(), "B", vapor,
            dry.gas().state().state().temperature, 400.0);
        let wet = dry.gas().exchange_species_boundary_relative(dry.catalog(),
            &[replacement], controls(), || false).unwrap();
        assert_eq!(wet.gas().model().thermodynamics().components(), model.thermodynamics().components());
        assert_eq!(wet.gas().model().water_species(), Some(&id("B")));
        near(wet.gas().species_mass_kg(&id("A")).unwrap(), carrier, 1e-14);
        near(wet.gas().species_mass_kg(&id("B")).unwrap(), vapor, 1e-14);
        near(wet.gas().state().state().temperature, 400.0, 2e-8);
        near(wet.gas().state().state().pressure, 100000.0, 2e-5);
        near(dry.step().energy_residual_j, 0.0, 1e-7);
        near(wet.step().energy_residual_j, 0.0, 1e-7);
    }
}

#[test]
fn relative_invalid_sources_inventory_and_cancelled_candidates_publish_nothing() {
    let model = catalog(0.0);
    let gas = parcel(&model);
    let saved = gas.clone();
    let valid = relative_stream(&model, "B", 1e-6, 400.0, 500.0);
    let vapor = gas.species_mass_kg(&id("B")).unwrap();
    let carrier = gas.species_mass_kg(&id("A")).unwrap();
    let bad = [
        vec![SpeciesBoundaryRelativeTransfer { relative_energy_into_gas_j: f64::NAN, ..valid.clone() }],
        vec![SpeciesBoundaryRelativeTransfer { mass_into_gas_kg: f64::INFINITY, ..valid.clone() }],
        vec![SpeciesBoundaryRelativeTransfer { species: id("foreign"), ..valid.clone() }],
        vec![valid.clone(), valid.clone()],
        vec![relative_stream(&model, "B", -2.0 * vapor, 400.0, 400.0)],
        vec![relative_stream(&model, "A", -carrier, 400.0, 400.0),
             relative_stream(&model, "B", -vapor, 400.0, 400.0)],
        vec![SpeciesBoundaryRelativeTransfer { relative_energy_into_gas_j: 1e9, ..valid.clone() }],
        vec![SpeciesBoundaryRelativeTransfer { mass_into_gas_kg: f64::MIN_POSITIVE, ..valid.clone() }],
    ];
    for rows in bad {
        assert!(gas.exchange_species_boundary_relative(&model, &rows, controls(), || false).is_err());
        assert_eq!(gas, saved);
    }
    let wrong_source = catalog(1.0);
    assert!(gas.exchange_species_boundary_relative(&wrong_source, &[valid.clone()], controls(), || false).is_err());
    let missing_source = model.with_composition(
        Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap(),
    ).unwrap();
    assert!(gas.exchange_species_boundary_relative(&missing_source, &[valid.clone()], controls(), || false).is_err());
    for temperature in [199.0, 2001.0, 0.0, f64::NAN] {
        assert!(model.species_specific_enthalpy_change_j_kg(&id("B"), temperature, temperature).is_err());
    }
    assert!(model.species_specific_enthalpy_change_j_kg(&id("foreign"), 400.0, 400.0).is_err());
    let mut calls = 0;
    gas.exchange_species_boundary_relative(&model, &[valid.clone()], controls(), || {
        calls += 1; false
    }).unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert_eq!(gas.exchange_species_boundary_relative(&model, &[valid.clone()], controls(), || {
            count += 1; count == stop
        }), Err(GasEnergyError::Cancelled));
        assert_eq!(gas, saved);
    }
}

#[test]
fn opposed_boundary_streams_cannot_hide_unrepresented_species_mass() {
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species_with_energy_reference("A", 0.03125, 3.5, 0.0, 0.0),
                 species_with_energy_reference("B", 0.03125, 3.5, 0.0, 0.0)],
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        ).unwrap(), vec![transport("A", 1.7e-5, 0.025), transport("B", 1.7e-5, 0.025)], None,
    ).unwrap();
    let density = parcel(&model).state().state().density;
    let gas = ClosedGasParcel::new(model.clone(), point(400.0, 100000.0), 2.0 / density).unwrap();
    assert_eq!(gas.species_mass_kg(&id("A")).unwrap(), 1.0);
    assert_eq!(gas.species_mass_kg(&id("B")).unwrap(), 1.0);
    let dm = 0.75 * f64::EPSILON;
    let rows = [relative_stream(&model, "A", dm, 400.0, 400.0),
                relative_stream(&model, "B", -dm, 400.0, 400.0)];
    let loose = gas.exchange_species_boundary_relative(&model, &rows, controls(), || false).unwrap();
    assert_eq!(loose.step().mass_residual_kg, 0.0);
    assert_eq!(loose.gas().species_mass_kg(&id("A")).unwrap(), 1.0 + f64::EPSILON);
    assert_eq!(loose.gas().species_mass_kg(&id("B")).unwrap(), 1.0 - f64::EPSILON);
    assert_eq!(loose.step().maximum_species_mass_residual_kg, f64::EPSILON / 4.0);
    let strict = GasSpeciesBoundarySolve { mass_tolerance_kg: f64::EPSILON / 8.0, ..controls() };
    let refusal = Err(GasEnergyError::Invalid("boundary gas endpoint did not retain individual species inventories"));
    assert_eq!(gas.exchange_species_boundary_relative(&model, &rows, strict, || false), refusal);
    let absolute = [stream(&model, "A", dm, 400.0), stream(&model, "B", -dm, 400.0)];
    assert_eq!(gas.exchange_species_boundary(&model, &absolute, strict, || false), refusal);
}
