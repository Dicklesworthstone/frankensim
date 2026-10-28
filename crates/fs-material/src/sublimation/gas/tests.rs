//! Consistent synthetic phase/caloric laws exercise the original gas owners.
use super::*;
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_qty::Temperature;
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

const MOLAR_MASS: f64 = 0.018015;
const LATENT: f64 = 200000.0;

#[derive(Clone, Copy)]
pub(super) struct Phase {
    pub(super) molar_mass: f64,
    pub(super) reference_temperature_k: f64,
}
impl SublimationLaw for Phase {
    fn temperature_range_k(&self) -> [f64; 2] {
        [200.0, 900.0]
    }
    fn properties(&self, temperature_k: f64) -> Result<SublimationProperties, SublimationError> {
        let specific_r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / self.molar_mass;
        let pressure = 1000.0
            * fs_math::det::exp(
                LATENT / specific_r * (1.0 / self.reference_temperature_k - 1.0 / temperature_k),
            );
        Ok(SublimationProperties {
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: pressure * LATENT
                / (specific_r * temperature_k * temperature_k),
            latent_heat_j_kg: LATENT,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 917.0,
            molar_mass_kg_mol: self.molar_mass,
        })
    }
}
pub(super) fn species() -> SpeciesId {
    SpeciesId::new("B").unwrap()
}
pub(super) fn catalog(reference: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.02897, 3.5, 0.0, 0.0),
                species_with_energy_reference("B", MOLAR_MASS, 4.0, 0.0, reference),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.99, 0.01]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.03)],
        Some(&species()),
    )
    .unwrap()
}
pub(super) fn gas(
    catalog: &ThermochemicalGasModel,
    temperature: f64,
    vapor_pressure: f64,
) -> ClosedGasParcel {
    let fraction = vapor_pressure / 100000.0;
    let model = catalog
        .with_composition(
            Composition::new(
                CompositionBasis::MoleFraction,
                vec![1.0 - fraction, fraction],
            )
            .unwrap(),
        )
        .unwrap();
    ClosedGasParcel::new(
        model,
        point(temperature, 100000.0).with("local-lot", 9.0).unwrap(),
        0.01,
    )
    .unwrap()
}
pub(super) fn binding(offset: f64) -> GasSurfaceBinding {
    GasSurfaceBinding {
        species: species(),
        accommodation_coefficient: 0.1,
        solid_enthalpy_offset_j_kg: offset,
    }
}
pub(super) fn controls() -> GasSurfaceSolve {
    GasSurfaceSolve {
        gas: GasSpeciesBoundarySolve {
            gas: crate::thermochemical_gas::energy::GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 900.0,
                energy_tolerance_j_kg: 1e-7,
                maximum_iterations: 128,
            },
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-7,
        },
        maximum_gross_inventory_fraction: 0.25,
        latent_compatibility_tolerance_j_kg: 1e-5,
    }
}

// The NASA-9 fixture is explicitly constant cp. Derive its caloric coefficients
// from the ORIGINAL species evaluation, then define h_s=h_v-L-offset exactly.
// This manufactured chart makes no empirical ice-accuracy claim.
pub(super) fn chart(catalog: &ThermochemicalGasModel, binding: &GasSurfaceBinding) -> (f64, f64) {
    let component = catalog
        .thermodynamics()
        .components()
        .iter()
        .find(|component| component.species() == &binding.species)
        .unwrap();
    let evaluated = component.evaluate(Temperature::new(400.0)).unwrap();
    let cp = evaluated.properties().cp().value() / component.molar_mass().value();
    let h = catalog
        .species_specific_enthalpy_j_kg(&binding.species, 400.0)
        .unwrap()
        - LATENT
        - binding.solid_enthalpy_offset_j_kg;
    (h, cp)
}
pub(super) fn temperature(h: f64, chart: (f64, f64)) -> Result<f64, SublimationError> {
    let t = 400.0 + (h - chart.0) / chart.1;
    if !t.is_finite() {
        Err(SublimationError::Invalid(
            "unrepresentable synthetic caloric inversion",
        ))
    } else {
        Ok(t)
    }
}
pub(super) fn solid(temperature_k: f64, chart: (f64, f64)) -> SolidSurfaceState {
    SolidSurfaceState {
        mass_kg: 0.01,
        area_m2: 1e-4,
        specific_enthalpy_j_kg: chart.0 + chart.1 * (temperature_k - 400.0),
    }
}
pub(super) fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}
pub(super) fn phase() -> Phase {
    Phase {
        molar_mass: MOLAR_MASS,
        reference_temperature_k: 400.0,
    }
}

#[test]
fn carrier_vacuum_sublimation_cools_solid_and_reconstructs_actual_gas_energy() {
    let catalog = catalog(0.0);
    let binding = binding(0.0);
    let chart = chart(&catalog, &binding);
    let solid = solid(400.0, chart);
    let gas = gas(&catalog, 400.0, 0.0);
    let carrier = gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap();
    let trial = advance_surface_with_gas(
        &phase(),
        &binding,
        solid,
        &gas,
        &catalog,
        |h| temperature(h, chart),
        1.0,
        controls(),
        || false,
    )
    .unwrap();
    assert!(trial.step().transferred_mass_kg > 0.0);
    assert_eq!(trial.step().incoming_mass_kg, 0.0);
    assert_eq!(trial.step().incoming_specific_energy_j_kg, None);
    assert!(trial.step().surface_temperature_k[1] < 400.0);
    assert!(trial.gas().species_partial_pressure_pa(&species()).unwrap() > 0.0);
    assert_ne!(
        trial.gas().state().state().temperature,
        gas.state().state().temperature
    );
    near(
        trial
            .gas()
            .species_mass_kg(&SpeciesId::new("A").unwrap())
            .unwrap(),
        carrier,
        1e-14,
    );
    near(trial.step().energy_residual_j, 0.0, 1e-7);
    near(trial.step().mass_residual_kg, 0.0, 1e-14);
    assert_eq!(trial.gas().point().axes().get("local-lot"), Some(&9.0));
    assert_eq!(trial.gas().model().water_species(), Some(&species()));
}

#[test]
fn deposition_heats_solid_and_both_gross_streams_enter_one_energy_account() {
    let catalog = catalog(-1.0e7);
    let binding = binding(-1.0e7);
    let chart = chart(&catalog, &binding);
    let solid = solid(400.0, chart);
    let gas = gas(&catalog, 400.0, 5000.0);
    let trial = advance_surface_with_gas(
        &phase(),
        &binding,
        solid,
        &gas,
        &catalog,
        |h| temperature(h, chart),
        0.2,
        controls(),
        || false,
    )
    .unwrap();
    assert!(trial.step().transferred_mass_kg < 0.0);
    assert!(trial.step().outgoing_mass_kg > 0.0);
    assert!(trial.step().incoming_mass_kg > trial.step().outgoing_mass_kg);
    assert!(trial.solid().mass_kg > solid.mass_kg);
    assert!(trial.step().surface_temperature_k[1] > 400.0);
    assert!(trial.step().recession_m < 0.0);
    near(
        trial.step().energy_into_gas_j,
        trial.step().transferred_mass_kg * trial.step().outgoing_specific_energy_j_kg,
        1e-8,
    );
    near(trial.step().energy_residual_j, 0.0, 1e-7);
    near(trial.step().mass_residual_kg, 0.0, 1e-14);
}

#[test]
fn zero_net_mass_at_unequal_temperatures_still_transfers_molecular_heat() {
    let catalog = catalog(0.0);
    let binding = binding(0.0);
    let chart = chart(&catalog, &binding);
    let solid = solid(800.0, chart);
    // Anchor this equally valid constant-latent pressure chart at 800 K so
    // the manufactured zero net flux is exact in floating-point arithmetic.
    let law = Phase {
        reference_temperature_k: 800.0,
        ..phase()
    };
    let pressure = law.properties(800.0).unwrap().saturation_pressure_pa * 0.5;
    let gas = gas(&catalog, 200.0, pressure);
    let trial = advance_surface_with_gas(
        &law,
        &binding,
        solid,
        &gas,
        &catalog,
        |h| temperature(h, chart),
        0.1,
        controls(),
        || false,
    )
    .unwrap();
    near(trial.step().transferred_mass_kg, 0.0, 1e-18);
    assert!(trial.step().outgoing_mass_kg > 0.0);
    assert!(trial.step().incoming_mass_kg > 0.0);
    assert!(trial.step().energy_into_gas_j > 0.0);
    assert!(trial.step().surface_temperature_k[1] < 800.0);
    assert!(trial.step().gas_temperature_k[1] > 200.0);
    near(trial.step().energy_residual_j, 0.0, 1e-7);
    // A net-mass limiter alone would miss this invalid counterflow timestep.
    assert!(
        advance_surface_with_gas(
            &law,
            &binding,
            solid,
            &gas,
            &catalog,
            |h| temperature(h, chart),
            100.0,
            controls(),
            || false
        )
        .is_err()
    );
}

#[test]
fn consistent_formation_reference_shifts_preserve_both_temperatures_and_mass() {
    let base = catalog(0.0);
    let shifted = catalog(1.0e7);
    let catalogs = [&base, &shifted];
    let bindings = [binding(0.0), binding(1.0e7)];
    let charts = [chart(&base, &bindings[0]), chart(&shifted, &bindings[1])];
    let gases = [gas(&base, 350.0, 500.0), gas(&shifted, 350.0, 500.0)];
    let trials = (0..2)
        .map(|index| {
            advance_surface_with_gas(
                &phase(),
                &bindings[index],
                solid(400.0, charts[index]),
                &gases[index],
                catalogs[index],
                |h| temperature(h, charts[index]),
                0.1,
                controls(),
                || false,
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    near(trials[0].solid().mass_kg, trials[1].solid().mass_kg, 1e-14);
    near(
        trials[0].step().surface_temperature_k[1],
        trials[1].step().surface_temperature_k[1],
        1e-7,
    );
    near(
        trials[0].step().gas_temperature_k[1],
        trials[1].step().gas_temperature_k[1],
        1e-7,
    );
    near(
        trials[1].step().energy_into_gas_j - trials[0].step().energy_into_gas_j,
        1.0e7 * trials[0].step().transferred_mass_kg,
        1e-8,
    );
    for trial in trials {
        near(trial.step().energy_residual_j, 0.0, 1e-7);
    }
}

#[test]
fn molecular_latent_and_sensible_targets_remain_analytic_at_large_species_datums() {
    let base = catalog(0.0);
    let caloric = chart(&base, &binding(0.0));
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MOLAR_MASS;
    for vapor_pressure in [500.0, 5000.0] {
        let mut original: Option<(f64, f64, f64)> = None;
        for reference in [0.0, -1e14, 1e14] {
            let source = catalog(reference);
            let alignment = binding(reference);
            let solid = solid(400.0, caloric);
            let gas = gas(&source, 350.0, vapor_pressure);
            let mut solve = controls();
            // The explicit source-alignment budget admits the ULP of the
            // declared datum; caloric energy and mass budgets stay unchanged.
            solve.latent_compatibility_tolerance_j_kg = 0.125;
            let trial = advance_surface_with_gas(&phase(), &alignment, solid, &gas, &source,
                |h| temperature(h, caloric), 0.1, solve, || false).unwrap();
            let report = trial.step();
            let dm = report.transferred_mass_kg;
            let outgoing = report.outgoing_mass_kg;
            let incoming = report.incoming_mass_kg;
            let alignment_residual = source.species_specific_enthalpy_j_kg(&species(), 400.0).unwrap()
                - (solid.specific_enthalpy_j_kg + alignment.solid_enthalpy_offset_j_kg) - LATENT;
            let molecular_difference = 3.5 * r * (400.0 - 350.0);
            let solid_heat = -dm * (LATENT + alignment_residual - 0.5 * r * 400.0)
                - incoming * molecular_difference;
            let expected_solid = 400.0 + solid_heat / ((solid.mass_kg - dm) * 4.0 * r);
            let gas_relative = outgoing * 4.0 * r * (400.0 - 350.0)
                - 0.5 * r * (outgoing * 400.0 - incoming * 350.0);
            let carrier = gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap();
            let vapor = gas.species_mass_kg(&species()).unwrap() + dm;
            let gas_capacity = carrier * 2.5 * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897
                + vapor * 3.0 * r;
            let expected_gas = 350.0 + (dm * r * 350.0 + gas_relative) / gas_capacity;
            near(report.surface_temperature_k[1], expected_solid, 2e-8);
            near(report.gas_temperature_k[1], expected_gas, 2e-8);
            near(trial.gas().state().state().pressure,
                (carrier / 0.02897 + vapor / MOLAR_MASS) * UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K
                    * expected_gas / gas.volume_m3(), 2e-5);
            near(report.energy_residual_j, 0.0, 1e-7);
            near(report.mass_residual_kg, 0.0, 1e-14);
            let state = (report.surface_temperature_k[1], report.gas_temperature_k[1], dm);
            if let Some(first) = original {
                let admitted_source_precision = (dm * alignment_residual).abs()
                    / ((solid.mass_kg - dm) * 4.0 * r);
                near(state.0, first.0, admitted_source_precision + 1e-10);
                near(state.1, first.1, 2e-8);
                assert_eq!(state.2, first.2);
            } else {
                original = Some(state);
            }
        }
    }
}

#[test]
fn zero_net_molecular_counterflow_retains_thermal_energy_beside_large_datums() {
    let base = catalog(0.0);
    let caloric = chart(&base, &binding(0.0));
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MOLAR_MASS;
    let law = Phase { reference_temperature_k: 800.0, ..phase() };
    for reference in [0.0, -1e14, 1e14] {
        let source = catalog(reference);
        let alignment = binding(reference);
        let solid = solid(800.0, caloric);
        let gas = gas(&source, 200.0, 500.0);
        let mut solve = controls();
        solve.latent_compatibility_tolerance_j_kg = 0.125;
        let trial = advance_surface_with_gas(&law, &alignment, solid, &gas, &source,
            |h| temperature(h, caloric), 0.1, solve, || false).unwrap();
        let report = trial.step();
        assert_eq!(report.transferred_mass_kg, 0.0);
        assert!(report.outgoing_mass_kg > 0.0 && report.incoming_mass_kg > 0.0);
        let heat = report.incoming_mass_kg * 3.5 * r * (800.0 - 200.0);
        near(report.energy_into_gas_j, heat, 1e-12);
        near(report.surface_temperature_k[1], 800.0 - heat / (solid.mass_kg * 4.0 * r), 2e-8);
        near(report.energy_residual_j, 0.0, 1e-7);
        assert_eq!(trial.gas().mass_kg(), gas.mass_kg());
        assert_eq!(trial.solid().mass_kg, solid.mass_kg);
        assert!(report.gas_temperature_k[1] > 200.0);
    }
}

#[test]
fn admitted_latent_alignment_error_remains_in_the_actual_solid_energy_target() {
    let source = catalog(0.0);
    let caloric = chart(&source, &binding(0.0));
    let alignment = binding(1.0);
    let solid = solid(400.0, caloric);
    let gas = gas(&source, 400.0, 0.0);
    let mut solve = controls();
    solve.latent_compatibility_tolerance_j_kg = 2.0;
    let trial = advance_surface_with_gas(&phase(), &alignment, solid, &gas, &source,
        |h| temperature(h, caloric), 1.0, solve, || false).unwrap();
    let delta = source.species_specific_enthalpy_j_kg(&species(), 400.0).unwrap()
        - (solid.specific_enthalpy_j_kg + alignment.solid_enthalpy_offset_j_kg) - LATENT;
    near(delta, -1.0, 1e-8);
    let dm = trial.step().transferred_mass_kg;
    assert!(dm * delta.abs() > 10.0 * solve.gas.energy_tolerance_j);
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MOLAR_MASS;
    let expected = 400.0 - dm * (LATENT + delta - 0.5 * r * 400.0)
        / ((solid.mass_kg - dm) * caloric.1);
    near(trial.step().surface_temperature_k[1], expected, 1e-10);
    // With this small datum, independent direct endpoint bookkeeping also
    // resolves the admitted alignment defect and must close the same first law.
    let actual_solid = trial.solid().mass_kg
        * (trial.solid().specific_enthalpy_j_kg + alignment.solid_enthalpy_offset_j_kg)
        - solid.mass_kg * (solid.specific_enthalpy_j_kg + alignment.solid_enthalpy_offset_j_kg);
    let actual_gas = trial.gas().mass_kg() * trial.gas().state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    near(actual_solid + actual_gas, 0.0, solve.gas.energy_tolerance_j);
    near(trial.step().energy_residual_j, 0.0, solve.gas.energy_tolerance_j);
}

#[test]
fn physical_timestep_refinement_converges_with_evolving_backpressure() {
    let catalog = catalog(0.0);
    let binding = binding(0.0);
    let chart = chart(&catalog, &binding);
    let advance = |steps: usize| {
        let mut solid = solid(400.0, chart);
        let mut gas = gas(&catalog, 400.0, 0.0);
        for _ in 0..steps {
            let candidate = advance_surface_with_gas(
                &phase(),
                &binding,
                solid,
                &gas,
                &catalog,
                |h| temperature(h, chart),
                2.0 / steps as f64,
                controls(),
                || false,
            )
            .unwrap();
            let (next_solid, next_gas, _, report) = candidate.into_parts();
            near(report.energy_residual_j, 0.0, 1e-7);
            solid = next_solid;
            gas = next_gas;
        }
        gas.species_mass_kg(&species()).unwrap()
    };
    let reference = advance(128);
    let errors = [1, 2, 4].map(|steps| (advance(steps) - reference).abs());
    assert!(
        errors[2] < 0.6 * errors[1] && errors[1] < 0.6 * errors[0],
        "{errors:?}"
    );
}

#[test]
fn inconsistent_calorics_molar_mass_budget_and_cancellation_preserve_both_inputs() {
    let catalog = catalog(0.0);
    let binding = binding(0.0);
    let chart = chart(&catalog, &binding);
    let solid = solid(400.0, chart);
    let gas = gas(&catalog, 400.0, 5000.0);
    let saved = gas.clone();
    let wrong_datum = GasSurfaceBinding {
        solid_enthalpy_offset_j_kg: 1000.0,
        ..binding.clone()
    };
    assert!(matches!(
        advance_surface_with_gas(
            &phase(),
            &wrong_datum,
            solid,
            &gas,
            &catalog,
            |h| temperature(h, chart),
            0.1,
            controls(),
            || false
        ),
        Err(GasSurfaceError::CaloricMismatch { .. })
    ));
    assert!(matches!(
        advance_surface_with_gas(
            &phase(),
            &binding,
            solid,
            &gas,
            &catalog,
            |h| temperature(h, chart).map(|t| if h == solid.specific_enthalpy_j_kg {
                t
            } else {
                t + 0.5
            }),
            0.1,
            controls(),
            || false
        ),
        Err(GasSurfaceError::CaloricMismatch { .. })
    ));
    assert!(
        advance_surface_with_gas(
            &Phase {
                molar_mass: 0.02,
                ..phase()
            },
            &binding,
            solid,
            &gas,
            &catalog,
            |h| temperature(h, chart),
            0.1,
            controls(),
            || false
        )
        .is_err()
    );
    assert!(
        advance_surface_with_gas(
            &phase(),
            &binding,
            solid,
            &gas,
            &catalog,
            |h| temperature(h, chart),
            10.0,
            controls(),
            || false
        )
        .is_err()
    );
    let mut restricted = controls();
    restricted.gas.gas.maximum_temperature_k = 300.0;
    assert!(
        advance_surface_with_gas(
            &phase(),
            &binding,
            solid,
            &gas,
            &catalog,
            |h| temperature(h, chart),
            0.1,
            restricted,
            || false
        )
        .is_err()
    );
    let mut calls = 0;
    advance_surface_with_gas(
        &phase(),
        &binding,
        solid,
        &gas,
        &catalog,
        |h| temperature(h, chart),
        0.1,
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        let result = advance_surface_with_gas(
            &phase(),
            &binding,
            solid,
            &gas,
            &catalog,
            |h| temperature(h, chart),
            0.1,
            controls(),
            || {
                count += 1;
                count == stop
            },
        );
        assert!(matches!(
            result,
            Err(GasSurfaceError::Surface(SublimationError::Cancelled))
                | Err(GasSurfaceError::Gas(GasEnergyError::Cancelled))
        ));
        assert_eq!(gas, saved);
    }
    assert_eq!(solid, self::solid(400.0, chart));
}
