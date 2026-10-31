//! Manufactured phase/caloric laws prove the coupled account independently.
mod water;

use super::*;
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_thermal::evaporation::EvaporationProperties;
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

const M: f64 = 0.018015;
const R: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / M;
const LATENT: f64 = 200000.0;
const SLOPE: f64 = 0.001;

struct Phase;
impl EvaporationLaw for Phase {
    fn temperature_range_k(&self) -> [f64; 2] {
        [200.0, 900.0]
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        let p = 1500.0 * fs_math::det::exp(LATENT / R * (1.0 / 300.0 - 1.0 / t));
        Ok(EvaporationProperties {
            saturation_pressure_pa: p,
            saturation_pressure_derivative_pa_k: p * LATENT / (R * t * t),
            latent_heat_j_kg: LATENT,
            liquid_density_kg_m3: 1000.0,
            molar_mass_kg_mol: M,
        })
    }
}
fn id() -> SpeciesId {
    SpeciesId::new("B").unwrap()
}
fn carrier() -> SpeciesId {
    SpeciesId::new("A").unwrap()
}
fn catalog(reference: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.02897, 3.5, 0.0, 0.0),
                species_with_energy_reference("B", M, 4.0, SLOPE, reference),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.99, 0.01]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.03)],
        Some(&id()),
    )
    .unwrap()
}
fn gas(catalog: &ThermochemicalGasModel, t: f64, pv: f64, volume: f64) -> ClosedGasParcel {
    let x = pv / 100000.0;
    ClosedGasParcel::new(
        catalog
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, vec![1.0 - x, x]).unwrap(),
            )
            .unwrap(),
        point(t, 100000.0).with("local-lot", 11.0).unwrap(),
        volume,
    )
    .unwrap()
}
fn liquid(t: f64) -> LiquidSurfaceState {
    LiquidSurfaceState {
        mass_kg: 0.01,
        area_m2: 1e-3,
        specific_enthalpy_j_kg: R * (4.0 * t + 0.5 * SLOPE * t * t) - LATENT,
    }
}
fn temperature(h: f64) -> Result<f64, EvaporationError> {
    let c = (h + LATENT) / R;
    Ok(2.0 * c / (4.0 + (16.0 + 2.0 * SLOPE * c).sqrt()))
}
fn binding(reference: f64) -> GasLiquidBinding {
    GasLiquidBinding {
        species: id(),
        liquid_enthalpy_offset_j_kg: reference,
    }
}
fn controls() -> GasEvaporationSolve {
    GasEvaporationSolve {
        gas: GasSpeciesBoundarySolve {
            gas: crate::thermochemical_gas::energy::GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 900.0,
                energy_tolerance_j_kg: 1e-8,
                maximum_iterations: 128,
            },
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-8,
        },
        maximum_inventory_fraction: 0.25,
        maximum_temperature_increment_k: 10.0,
        latent_compatibility_tolerance_j_kg: 1e-5,
    }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a} != {b}; tolerance {tol}");
}

#[test]
fn dry_carrier_receives_vapor_and_actual_transport_pressure_and_temperature_change() {
    let catalog = catalog(0.0);
    let gas = gas(&catalog, 300.0, 0.0, 0.01);
    let liquid = liquid(300.0);
    let before = gas.clone();
    let step = advance_surface_with_gas(
        &Phase,
        &binding(0.0),
        liquid,
        &gas,
        &catalog,
        temperature,
        0.02,
        1.0,
        controls(),
        || false,
    )
    .unwrap();
    let report = step.step();
    assert!(report.transferred_mass_kg > 0.0);
    assert!(step.liquid().mass_kg < liquid.mass_kg);
    assert!(report.surface_temperature_k[1] < 300.0);
    assert!(report.vapor_partial_pressure_pa[1] > 0.0);
    assert_ne!(
        step.gas().state().state().pressure,
        gas.state().state().pressure
    );
    assert_ne!(
        step.gas().state().state().temperature,
        gas.state().state().temperature
    );
    assert_ne!(step.gas().state(), gas.state());
    near(
        step.gas().species_mass_kg(&carrier()).unwrap(),
        gas.species_mass_kg(&carrier()).unwrap(),
        1e-14,
    );
    near(
        step.gas().species_mass_kg(&id()).unwrap(),
        report.transferred_mass_kg,
        1e-14,
    );
    near(
        report.liquid_enthalpy_change_j + report.gas_internal_energy_change_j,
        0.0,
        1e-8,
    );
    near(report.energy_residual_j, 0.0, 1e-8);
    near(report.mass_residual_kg, 0.0, 1e-14);
    assert_eq!(step.gas().point().axes().get("local-lot"), Some(&11.0));
    assert_eq!(step.gas().model().water_species(), Some(&id()));
    assert_eq!(step.gas().volume_m3(), gas.volume_m3());
    assert!(std::ptr::eq(step.catalog(), &catalog));
    assert_eq!(gas, before);
}

#[test]
fn cold_surface_condenses_upstream_vapor_with_latent_and_sensible_energy_once() {
    let catalog = catalog(-1e7);
    let gas = gas(&catalog, 320.0, 5000.0, 0.01);
    let liquid = liquid(300.0);
    let step = advance_surface_with_gas(
        &Phase,
        &binding(-1e7),
        liquid,
        &gas,
        &catalog,
        temperature,
        0.02,
        0.5,
        controls(),
        || false,
    )
    .unwrap();
    let report = step.step();
    assert!(report.transferred_mass_kg < 0.0);
    assert!(step.liquid().mass_kg > liquid.mass_kg);
    assert!(report.surface_temperature_k[1] > 300.0);
    assert!(report.vapor_partial_pressure_pa[1] < 5000.0);
    let upstream = catalog
        .species_specific_enthalpy_j_kg(&id(), 320.0)
        .unwrap();
    near(
        report.energy_into_gas_j,
        report.transferred_mass_kg * upstream,
        1e-8,
    );
    near(
        report.latent_energy_j,
        report.transferred_mass_kg * LATENT,
        1e-12,
    );
    let expected = -report.transferred_mass_kg
        * (LATENT + R * (4.0 * 20.0 + 0.5 * SLOPE * (320.0 * 320.0 - 300.0 * 300.0)));
    let sensible = step.liquid().mass_kg
        * (step.liquid().specific_enthalpy_j_kg - liquid.specific_enthalpy_j_kg);
    near(sensible, expected, 1e-8);
    near(
        report.liquid_enthalpy_change_j + report.gas_internal_energy_change_j,
        0.0,
        1e-8,
    );
}

#[test]
fn finite_backpressure_reduces_evaporation_and_frozen_rate_overshoot_refuses() {
    let catalog = catalog(0.0);
    let mut current = gas(&catalog, 300.0, 0.0, 1e-4);
    let mut surface = liquid(300.0);
    let mut previous = f64::INFINITY;
    for _ in 0..20 {
        let result = advance_surface_with_gas(
            &Phase,
            &binding(0.0),
            surface,
            &current,
            &catalog,
            temperature,
            0.02,
            0.05,
            controls(),
            || false,
        )
        .unwrap();
        assert!(result.step().transferred_mass_kg > 0.0);
        assert!(result.step().transferred_mass_kg < previous);
        previous = result.step().transferred_mass_kg;
        (surface, current, _, _) = result.into_parts();
    }
    let tiny = gas(&catalog, 300.0, 0.0, 1e-5);
    let original = tiny.clone();
    let mut config = controls();
    config.maximum_temperature_increment_k = 100.0;
    let error = advance_surface_with_gas(
        &Phase,
        &binding(0.0),
        liquid(300.0),
        &tiny,
        &catalog,
        temperature,
        0.02,
        2.0,
        config,
        || false,
    )
    .unwrap_err();
    assert!(
        matches!(error, GasEvaporationError::StepSizeLimit { .. }),
        "{error:?}"
    );
    assert_eq!(tiny, original);
}

#[test]
fn equilibrium_sealed_film_and_zero_duration_preserve_both_states() {
    let catalog = catalog(0.0);
    for (pv, coefficient, dt) in [(1500.0, 0.02, 1.0), (500.0, 0.0, 1.0), (500.0, 0.02, 0.0)] {
        let gas = gas(&catalog, 300.0, pv, 0.01);
        let liquid = liquid(300.0);
        let step = advance_surface_with_gas(
            &Phase,
            &binding(0.0),
            liquid,
            &gas,
            &catalog,
            temperature,
            coefficient,
            dt,
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(step.liquid(), liquid);
        assert_eq!(step.gas(), &gas);
        assert_eq!(step.step().transferred_mass_kg, 0.0);
        assert_eq!(step.step().energy_residual_j, 0.0);
    }
}

#[test]
fn datum_alignment_changes_reports_without_changing_physical_transfer() {
    let mut states = Vec::new();
    for reference in [0.0, -1e8, 1e12] {
        let catalog = catalog(reference);
        let gas = gas(&catalog, 310.0, 500.0, 0.01);
        let mut config = controls();
        config.latent_compatibility_tolerance_j_kg = 0.01;
        let step = advance_surface_with_gas(
            &Phase,
            &binding(reference),
            liquid(300.0),
            &gas,
            &catalog,
            temperature,
            0.02,
            0.5,
            config,
            || false,
        )
        .unwrap();
        states.push((
            step.liquid(),
            step.gas().state().state().temperature,
            step.step().transferred_mass_kg,
        ));
        near(step.step().energy_residual_j, 0.0, 1e-8);
    }
    for state in &states[1..] {
        near(state.0.mass_kg, states[0].0.mass_kg, 0.0);
        near(
            state.0.specific_enthalpy_j_kg,
            states[0].0.specific_enthalpy_j_kg,
            1e-8,
        );
        near(state.1, states[0].1, 1e-10);
        near(state.2, states[0].2, 0.0);
    }
}

#[test]
fn physical_time_refinement_converges_for_coupled_liquid_and_finite_atmosphere() {
    let catalog = catalog(0.0);
    let run = |n: usize| {
        let mut gas = gas(&catalog, 300.0, 200.0, 0.001);
        let mut liquid = liquid(300.0);
        for _ in 0..n {
            let result = advance_surface_with_gas(
                &Phase,
                &binding(0.0),
                liquid,
                &gas,
                &catalog,
                temperature,
                0.02,
                10.0 / n as f64,
                controls(),
                || false,
            )
            .unwrap();
            (liquid, gas, _, _) = result.into_parts();
        }
        gas.species_partial_pressure_pa(&id()).unwrap()
    };
    let reference = run(2048);
    let errors = [16, 32, 64].map(|n| (run(n) - reference).abs());
    assert!(errors[0] > 1e-3, "inactive refinement: {errors:?}");
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!(
            (0.8..1.2).contains(&order),
            "first-order ladder: {errors:?}, p={order}"
        );
    }
}

#[test]
fn bad_datum_species_inventory_and_late_cancellation_never_publish() {
    let catalog = catalog(0.0);
    let gas = gas(&catalog, 300.0, 5000.0, 0.01);
    let before = gas.clone();
    let surface = liquid(300.0);
    let run = |binding: &GasLiquidBinding, dt, config, cancelled: &mut dyn FnMut() -> bool| {
        advance_surface_with_gas(
            &Phase,
            binding,
            surface,
            &gas,
            &catalog,
            temperature,
            0.02,
            dt,
            config,
            cancelled,
        )
    };
    assert!(matches!(
        run(&binding(1e5), 1.0, controls(), &mut || false),
        Err(GasEvaporationError::CaloricMismatch { .. })
    ));
    let mut bad = binding(0.0);
    bad.species = carrier();
    assert!(matches!(
        run(&bad, 1.0, controls(), &mut || false),
        Err(GasEvaporationError::Invalid(_))
    ));
    let error = run(&binding(0.0), 10000.0, controls(), &mut || false).unwrap_err();
    assert!(matches!(error, GasEvaporationError::StepSizeLimit { .. }));
    let mut calls = 0;
    run(&binding(0.0), 1.0, controls(), &mut || {
        calls += 1;
        false
    })
    .unwrap();
    for stop in [1, 2, calls] {
        let mut count = 0;
        assert!(
            run(&binding(0.0), 1.0, controls(), &mut || {
                count += 1;
                count >= stop
            })
            .is_err()
        );
    }
    assert_eq!(gas, before);
}
