//! G0/G1 source feedback, original coupled conservation and transactional refusal.
//! Reuse the existing molecular-film NASA9/material fixture unchanged.
use super::*;
use crate::graph::thermal::medium::mass_transfer::evaporation::{
    EvaporationFilmError, SourceLaminarPlateEvaporationFilms,
};
use fs_material::evaporation::gas::implicit::adaptive::AdaptiveImplicitGasControls;
use fs_material::evaporation::gas::implicit::{ImplicitGasControls, ImplicitGasSurface};
use fs_material::evaporation::gas::{GasEvaporationSolve, GasLiquidBinding};
use fs_thermal::enthalpy::backward_euler::{
    BoundaryExchangeStepCoefficients, BoundaryStepCoefficients, StepCoefficients,
    TransportStepCoefficients,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermal::evaporation::{EvaporationError, EvaporationLaw, EvaporationProperties, GasFilm};

const LIQUID_MASS: f64 = 1e-4;
const LIQUID_AREA: f64 = 1e-3;

struct LiquidPhase;
impl EvaporationLaw for LiquidPhase {
    fn temperature_range_k(&self) -> [f64; 2] {
        [200.0, 900.0]
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        let source = Phase.properties(t).unwrap();
        Ok(EvaporationProperties {
            saturation_pressure_pa: source.saturation_pressure_pa,
            saturation_pressure_derivative_pa_k: source.saturation_pressure_derivative_pa_k,
            latent_heat_j_kg: LATENT,
            liquid_density_kg_m3: 1000.0,
            molar_mass_kg_mol: MOLAR_MASS,
        })
    }
}

struct LiquidCaloric {
    h400: f64,
    cp: f64,
}
impl LiquidCaloric {
    fn new(catalog: &ThermochemicalGasModel) -> Self {
        let (h400, cp) = chart(catalog);
        Self { h400, cp }
    }
    fn h(&self, t: f64) -> f64 {
        self.h400 + self.cp * (t - 400.0)
    }
}
impl ThermalLaw for LiquidCaloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 400.0 + (h - self.h400) / self.cp;
        if t.is_finite() && (200.0..=900.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp
    }
}

fn heat_graph<'a>(
    law: &'a LiquidCaloric,
    t: f64,
    power_w: f64,
    heater_off_s: Option<f64>,
) -> DepositionSimulation<'a> {
    let boundary = Boundary {
        ambient_k: 400.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: power_w,
    };
    let mut events = vec![TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id: CellId(4),
            segment_id: 4,
            body_id: 4,
            material_id: MaterialId(1),
            mass_kg: LIQUID_MASS,
            specific_h_j_kg: law.h(t),
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [0.1, 0.0, 0.0],
                volume_m3: LIQUID_MASS / 1000.0,
                surface_area_m2: LIQUID_AREA,
            },
            boundary,
        }),
    }];
    if let Some(time_s) = heater_off_s {
        events.push(TimedEvent {
            time_s,
            event: Event::Boundary(
                CellId(4),
                Boundary {
                    source_w: 0.0,
                    ..boundary
                },
            ),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 4096,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-6,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn liquid_binding() -> GasLiquidBinding {
    GasLiquidBinding {
        species: SpeciesId::new("B").unwrap(),
        liquid_enthalpy_offset_j_kg: 0.0,
    }
}
fn liquid_surface<'a>(
    law: &LiquidCaloric,
    binding: &'a GasLiquidBinding,
) -> ImplicitGasSurface<'a> {
    ImplicitGasSurface {
        cell_id: CellId(4),
        law: &LiquidPhase,
        binding,
        pure_liquid_enthalpy_range_j_kg: [law.h(250.0), law.h(500.0)],
        // A sealed authored default makes omission of live source feedback
        // fail the physical evaporation assertions, rather than pass silently.
        mass_transfer_coefficient_m_s: 0.0,
        heat_transfer_coefficient_w_m2_k: 0.0,
    }
}
fn phase_controls() -> ImplicitGasControls {
    let mut solve = ImplicitGasControls::new(GasEvaporationSolve {
        gas: GasSpeciesBoundarySolve {
            gas: GasEnergySolve {
                minimum_temperature_k: 250.0,
                maximum_temperature_k: 500.0,
                energy_tolerance_j_kg: 1e-8,
                maximum_iterations: 128,
            },
            mass_tolerance_kg: 1e-12,
            energy_tolerance_j: 1e-6,
        },
        maximum_inventory_fraction: 0.25,
        maximum_temperature_increment_k: 100.0,
        latent_compatibility_tolerance_j_kg: 1e-5,
    });
    solve.thermal.difference_step_j = 1e-4;
    solve
}
fn heat<T>(payload: T) -> BoundaryExchangeStepCoefficients<T> {
    TransportStepCoefficients::from(BoundaryStepCoefficients::from(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload,
    }))
    .into()
}
fn films(gas: &ClosedGasParcel) -> BTreeMap<CellId, GasFilm> {
    BTreeMap::from([(
        CellId(4),
        GasFilm {
            total_pressure_pa: gas.state().state().pressure,
            vapor_partial_pressure_pa: gas
                .species_partial_pressure_pa(&SpeciesId::new("B").unwrap())
                .unwrap(),
            temperature_k: f64::NAN,
            mass_transfer_coefficient_m_s: f64::NAN,
        },
    )])
}

#[test]
fn liquid_film_refreshes_original_diffusivity_viscosity_composition_and_temperature() {
    let catalog = catalog();
    let initial_gas = gas(&catalog, 400.0, 100000.0, 0.003);
    let changed_gas = gas(&catalog, 450.0, 130000.0, 0.04);
    let caloric = LiquidCaloric::new(&catalog);
    let thermal = heat_graph(&caloric, 350.0, 0.0, None);
    let liquid_binding = liquid_binding();
    let surfaces = [liquid_surface(&caloric, &liquid_binding)];
    let card = material(vec![scalar(
        BINARY_DIFFUSIVITY_PROPERTY,
        BINARY_DIFFUSIVITY_DIMS,
        "pressure",
        vec![(50000.0, 2e-5), (150000.0, 2e-5 / 3.0)],
        [200.0, 900.0],
    )]);
    let context = point(250.0, 50000.0);
    for convention in [
        ChamberFilmTemperature::Bulk,
        ChamberFilmTemperature::MeanSurfaceAndBulk,
    ] {
        let binding = LaminarPlateFilmBinding {
            property_temperature: convention,
            ..source(&card, &context)
        };
        let prepared = SourceLaminarPlateEvaporationFilms::new(
            &BTreeMap::from([(CellId(4), binding.clone())]),
            &surfaces,
            &thermal,
            &initial_gas,
            4,
            || false,
        )
        .unwrap();
        let initial_selection = prepared.selected_diffusivity(CellId(4)).unwrap();
        assert_eq!(
            initial_selection.state_coordinate("pressure"),
            Some(100000.0)
        );
        let expected_initial = match convention {
            ChamberFilmTemperature::Bulk => 400.0,
            ChamberFilmTemperature::MeanSurfaceAndBulk => 375.0,
        };
        assert_eq!(
            initial_selection.state_coordinate("temperature"),
            Some(expected_initial)
        );
        let mut original = films(&initial_gas);
        prepared
            .resolve(&thermal, &initial_gas, &mut original, || false)
            .unwrap();
        let original_coefficient = original[&CellId(4)].mass_transfer_coefficient_m_s;
        let mut current = films(&changed_gas);
        let total = current[&CellId(4)].total_pressure_pa.to_bits();
        let partial = current[&CellId(4)].vapor_partial_pressure_pa.to_bits();
        prepared
            .resolve(&thermal, &changed_gas, &mut current, || false)
            .unwrap();
        let film = current[&CellId(4)];
        assert_eq!(film.total_pressure_pa.to_bits(), total);
        assert_eq!(film.vapor_partial_pressure_pa.to_bits(), partial);
        let temperature = match convention {
            ChamberFilmTemperature::Bulk => 450.0,
            ChamberFilmTemperature::MeanSurfaceAndBulk => 400.0,
        };
        assert_eq!(film.temperature_k, temperature);
        let actual_fluid = changed_gas
            .model()
            .evaluate(&point(temperature, 130000.0))
            .unwrap();
        let expected = evaluate_laminar_plate_gas_film(
            &actual_fluid,
            &binding.source,
            &liquid_binding.species,
            LaminarPlateFilmConditions {
                start: Length::new(0.02),
                end: Length::new(0.1),
                speed_m_s: 1.0,
                temperature: Temperature::new(temperature),
                total_pressure: Pressure::new(130000.0),
            },
            || false,
        )
        .unwrap();
        near(
            film.mass_transfer_coefficient_m_s,
            expected.mass_transfer_coefficient_m_s(),
            1e-15,
        );
        assert_ne!(film.mass_transfer_coefficient_m_s, original_coefficient);
        // Speculative source evaluations do not rewrite the prepared selection.
        assert_eq!(
            prepared
                .selected_diffusivity(CellId(4))
                .unwrap()
                .state_coordinate("pressure"),
            Some(100000.0)
        );
        assert_eq!(binding.source.context, &context);
    }
}

#[test]
fn original_endpoint_evaporates_condenses_and_closes_actual_gas_liquid_energy() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let card = diffusion();
    let context = point(400.0, 100000.0);
    let binding = liquid_binding();
    let surfaces = [liquid_surface(&caloric, &binding)];
    for vapor_fraction in [0.003, 0.03] {
        let thermal = heat_graph(&caloric, 400.0, 0.0, None);
        let initial_gas = gas(&catalog, 400.0, 100000.0, vapor_fraction);
        let prepared = SourceLaminarPlateEvaporationFilms::new(
            &BTreeMap::from([(CellId(4), source(&card, &context))]),
            &surfaces,
            &thermal,
            &initial_gas,
            4,
            || false,
        )
        .unwrap();
        let duration = 1.0;
        let candidate = prepared
            .solve(
                &thermal,
                &initial_gas,
                &catalog,
                duration,
                phase_controls(),
                |_, _, _| Ok(heat(())),
                || false,
            )
            .unwrap();
        let liquid = &candidate.thermal().cells()[&CellId(4)];
        let dm = LIQUID_MASS - liquid.mass_kg;
        assert_eq!(dm > 0.0, vapor_fraction < 0.01);
        assert_eq!(
            candidate.thermal().temperature_k(CellId(4)).unwrap() < 400.0,
            dm > 0.0
        );
        near(
            candidate.gas().species_mass_kg(&binding.species).unwrap()
                - initial_gas.species_mass_kg(&binding.species).unwrap(),
            dm,
            1e-12,
        );
        near(
            candidate
                .gas()
                .species_mass_kg(&SpeciesId::new("A").unwrap())
                .unwrap(),
            initial_gas
                .species_mass_kg(&SpeciesId::new("A").unwrap())
                .unwrap(),
            1e-14,
        );
        let original_energy = LIQUID_MASS * caloric.h(400.0)
            + initial_gas.mass_kg() * initial_gas.state().specific_internal_energy_j_kg();
        let final_energy = liquid.mass_kg * liquid.specific_h_j_kg
            + candidate.gas().mass_kg() * candidate.gas().state().specific_internal_energy_j_kg();
        near(final_energy, original_energy, 1e-6);
        // Independently evaluate the original endpoint's full Stefan formula.
        let mut endpoint_films = films(candidate.gas());
        prepared
            .resolve(
                candidate.thermal(),
                candidate.gas(),
                &mut endpoint_films,
                || false,
            )
            .unwrap();
        let f = endpoint_films[&CellId(4)];
        let tl = candidate.thermal().temperature_k(CellId(4)).unwrap();
        let ps = 1000.0 * (LATENT / (R / MOLAR_MASS) * (1.0 / 400.0 - 1.0 / tl)).exp();
        let expected = duration
            * LIQUID_AREA
            * MOLAR_MASS
            * f.total_pressure_pa
            * f.mass_transfer_coefficient_m_s
            / (R * f.temperature_k)
            * ((f.total_pressure_pa - f.vapor_partial_pressure_pa) / (f.total_pressure_pa - ps))
                .ln();
        near(dm, expected, 1e-12);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells()[&CellId(4)].mass_kg, LIQUID_MASS);
    }
}

#[test]
fn source_velocity_changes_real_evaporation_instead_of_an_unused_coefficient() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let thermal = heat_graph(&caloric, 400.0, 0.0, None);
    let initial_gas = gas(&catalog, 400.0, 100000.0, 0.0);
    let card = diffusion();
    let context = point(400.0, 100000.0);
    let binding = liquid_binding();
    let surfaces = [liquid_surface(&caloric, &binding)];
    let mut masses = Vec::new();
    for speed in [1.0, 4.0] {
        let mut film = source(&card, &context);
        film.geometry.speed_m_s = speed;
        let prepared = SourceLaminarPlateEvaporationFilms::new(
            &BTreeMap::from([(CellId(4), film)]),
            &surfaces,
            &thermal,
            &initial_gas,
            4,
            || false,
        )
        .unwrap();
        let candidate = prepared
            .solve(
                &thermal,
                &initial_gas,
                &catalog,
                0.1,
                phase_controls(),
                |_, _, _| Ok(heat(())),
                || false,
            )
            .unwrap();
        let dm = LIQUID_MASS - candidate.thermal().cells()[&CellId(4)].mass_kg;
        assert!(dm > 0.0);
        assert!(candidate.gas().species_mass_kg(&binding.species).unwrap() > 0.0);
        masses.push(dm);
    }
    assert!(masses[1] > 1.99 * masses[0] && masses[1] < 2.01 * masses[0]);
    assert_eq!(initial_gas.species_mass_kg(&binding.species).unwrap(), 0.0);
}

#[test]
fn adaptive_source_feedback_respects_heater_event_and_whole_window_conservation() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let thermal = heat_graph(&caloric, 400.0, 0.1, Some(1.0));
    let initial_gas = gas(&catalog, 400.0, 100000.0, 0.003);
    let card = diffusion();
    let context = point(400.0, 100000.0);
    let binding = liquid_binding();
    let surfaces = [liquid_surface(&caloric, &binding)];
    let prepared = SourceLaminarPlateEvaporationFilms::new(
        &BTreeMap::from([(CellId(4), source(&card, &context))]),
        &surfaces,
        &thermal,
        &initial_gas,
        4,
        || false,
    )
    .unwrap();
    let mut controls = AdaptiveImplicitGasControls::new(phase_controls());
    controls.initial_dt_s = 0.5;
    controls.maximum_dt_s = 0.5;
    controls.solve.thermal.maximum_evaluations = 100000;
    controls.solve.thermal.maximum_iterations = 10000;
    controls.temperature_absolute_tolerance_k = 0.01;
    controls.enthalpy_absolute_tolerance_j_kg = 20.0;
    controls.mass_absolute_tolerance_kg = 1e-9;
    let result = prepared
        .advance(
            &thermal,
            &initial_gas,
            &catalog,
            2.0,
            controls,
            |thermal, _, _| Ok(heat(thermal.time_s())),
            || false,
        )
        .unwrap();
    assert_eq!(result.thermal().time_s(), 2.0);
    assert_eq!(result.payload(), Some(&2.0));
    assert_eq!(result.thermal().cells()[&CellId(4)].boundary.source_w, 0.0);
    near(result.report().external_heat_j, 0.1, 1e-8);
    let final_liquid = &result.thermal().cells()[&CellId(4)];
    assert!(final_liquid.mass_kg < LIQUID_MASS);
    let energy_change = final_liquid.mass_kg * final_liquid.specific_h_j_kg
        - LIQUID_MASS * caloric.h(400.0)
        + result.gas().mass_kg() * result.gas().state().specific_internal_energy_j_kg()
        - initial_gas.mass_kg() * initial_gas.state().specific_internal_energy_j_kg();
    near(energy_change, 0.1, 1e-6);
    near(
        final_liquid.mass_kg + result.gas().mass_kg(),
        LIQUID_MASS + initial_gas.mass_kg(),
        1e-12,
    );
    assert!(
        result.report().absolute_energy_residual_j
            <= controls.solve.thermal.total_energy_tolerance_j
    );
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn liquid_source_exhaustion_cannot_switch_claims_or_publish_partial_endpoint() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let thermal = heat_graph(&caloric, 400.0, 1.0, None);
    let original_cells = thermal.cells().clone();
    let initial_gas = gas(&catalog, 400.0, 100000.0, 0.003);
    let original_gas = initial_gas.clone();
    let card = material(vec![
        scalar(
            BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS,
            "temperature",
            vec![(200.0, 1e-5), (400.01, 1e-5)],
            [200.0, 400.01],
        ),
        scalar(
            BINARY_DIFFUSIVITY_PROPERTY,
            BINARY_DIFFUSIVITY_DIMS,
            "temperature",
            vec![(400.02, 2e-5), (900.0, 2e-5)],
            [400.02, 900.0],
        ),
    ]);
    let context = point(400.0, 100000.0);
    let binding = liquid_binding();
    let surfaces = [liquid_surface(&caloric, &binding)];
    let prepared = SourceLaminarPlateEvaporationFilms::new(
        &BTreeMap::from([(CellId(4), source(&card, &context))]),
        &surfaces,
        &thermal,
        &initial_gas,
        4,
        || false,
    )
    .unwrap();
    let result = prepared.solve(
        &thermal,
        &initial_gas,
        &catalog,
        1.0,
        phase_controls(),
        |_, _, _| Ok(heat(())),
        || false,
    );
    assert!(matches!(
        result,
        Err(EvaporationFilmError::Film(GasFilmError::Source(_)))
    ));
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(initial_gas, original_gas);
    let cancelled = prepared.solve(
        &thermal,
        &initial_gas,
        &catalog,
        0.1,
        phase_controls(),
        |_, _, _| Ok(heat(())),
        || true,
    );
    assert!(cancelled.is_err());
    assert_eq!(thermal.cells(), &original_cells);
    assert_eq!(initial_gas, original_gas);
}

#[test]
fn liquid_plate_domain_excludes_strong_loading_foreign_pairs_and_stagnant_flow() {
    let catalog = catalog();
    let caloric = LiquidCaloric::new(&catalog);
    let thermal = heat_graph(&caloric, 400.0, 0.0, None);
    let initial_gas = gas(&catalog, 400.0, 100000.0, 0.003);
    let loaded = gas(&catalog, 400.0, 100000.0, 0.06);
    let card = diffusion();
    let context = point(400.0, 100000.0);
    let binding = liquid_binding();
    let surfaces = [liquid_surface(&caloric, &binding)];
    let bindings = BTreeMap::from([(CellId(4), source(&card, &context))]);
    assert!(matches!(
        SourceLaminarPlateEvaporationFilms::new(&bindings, &surfaces, &thermal, &loaded, 4, || {
            false
        },),
        Err(EvaporationFilmError::OutsideDilute { .. })
    ));
    let mut invalid_pair = source(&card, &context);
    invalid_pair.source.species[0] = SpeciesId::new("C").unwrap();
    assert!(matches!(
        SourceLaminarPlateEvaporationFilms::new(
            &BTreeMap::from([(CellId(4), invalid_pair)]),
            &surfaces,
            &thermal,
            &initial_gas,
            4,
            || false,
        ),
        Err(EvaporationFilmError::Film(GasFilmError::Invalid(_)))
    ));
    let mut stagnant = source(&card, &context);
    stagnant.geometry.speed_m_s = 0.0;
    assert!(matches!(
        SourceLaminarPlateEvaporationFilms::new(
            &BTreeMap::from([(CellId(4), stagnant)]),
            &surfaces,
            &thermal,
            &initial_gas,
            4,
            || false,
        ),
        Err(EvaporationFilmError::Film(GasFilmError::Invalid(_)))
    ));
    assert!(matches!(
        SourceLaminarPlateEvaporationFilms::new(
            &bindings,
            &surfaces,
            &thermal,
            &initial_gas,
            4,
            || true,
        ),
        Err(EvaporationFilmError::Film(GasFilmError::Cancelled))
    ));
}
