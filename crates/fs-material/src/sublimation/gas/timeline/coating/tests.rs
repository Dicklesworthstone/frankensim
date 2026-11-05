//! Original thermal cells and NASA gas owners; manufactured compatible charts.
use super::*;
use crate::sublimation::gas::{chamber::GasChamberSolve, tests as fixture};
use fs_thermal::enthalpy::sublimation::coating::CoatingSurface;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, DepositionSimulation, MaterialId,
    ThermalLaw, ThermalLink, TimedEvent,
};
use fs_thermal::sublimation::{SolverConfig, SurfaceKinetics};
use fs_thermochem::SpeciesId;

#[derive(Clone, Copy)]
struct Substrate {
    datum: f64,
}
impl ThermalLaw for Substrate {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + (h - self.datum) / 800.0;
        if (200.0..=900.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 800.0
    }
}

fn cell(id: u64, chart: Substrate, power: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(1),
        mass_kg: 0.01,
        specific_h_j_kg: chart.datum,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 0.01, 0.0],
            volume_m3: 0.01 / 2700.0,
            surface_area_m2: 1e-4,
        },
        boundary: Boundary {
            ambient_k: 350.0,
            htc_w_m2_k: 0.0,
            source_w: power,
            reservoir: None,
        },
    }
}
fn birth(time_s: f64, cell: DepositedCell) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Birth(cell),
    }
}
fn controls() -> GasSupportedChamberSolve {
    GasSupportedChamberSolve {
        chamber: GasChamberSolve {
            surface: fixture::controls(),
            maximum_surfaces: 2,
        },
        storage: SolverConfig {
            relative_tolerance: 1e-12,
            absolute_energy_tolerance_j: 1e-9,
            max_iterations: 128,
        },
    }
}
fn heat(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-7,
    }
}
fn selected(bindings: &[GasSurfaceBinding; 2]) -> [GasCoatingTimelineSurface<'_>; 2] {
    std::array::from_fn(|i| GasCoatingTimelineSurface {
        cell_id: CellId(i as u64 + 1),
        binding: &bindings[i],
    })
}
fn initial<'a>(
    substrate: &'a Substrate,
    law: &'a dyn fs_thermal::sublimation::SublimationLaw,
    deposit_h: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
    bindings: &[GasSurfaceBinding; 2],
    events: Vec<TimedEvent>,
    masses: [f64; 2],
    max_steps: usize,
) -> CoatingSimulation<'a> {
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(1), substrate)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 8,
            max_substeps: max_steps,
        },
    )
    .unwrap();
    let surfaces: [_; 2] = std::array::from_fn(|i| CoatingSurface {
        cell_id: CellId(i as u64 + 1),
        law,
        deposit_enthalpy: deposit_h,
        kinetics: SurfaceKinetics {
            accommodation_coefficient: bindings[i].accommodation_coefficient,
            mass_transfer_coefficient_m_s: None,
        },
        initial_mass_kg: masses[i],
        reservoir: VaporReservoir::Prescribed {
            partial_pressure_pa: 6000.0,
            temperature_k: 350.0,
        },
    });
    CoatingSimulation::new(thermal, &surfaces, || Ok(())).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e} (tolerance {tolerance:e})"
    );
}
fn check(candidate: &GasCoatingTimelineCandidate<'_, '_>) {
    let audit = candidate.audit();
    close(audit.mass_residual_kg, 0.0, 1e-14);
    close(audit.energy_residual_j, 0.0, 1e-7);
    assert!(audit.condensed_phases.is_empty());
    let gas = candidate.gas();
    for (&id, coating) in candidate.coating().coatings() {
        assert_eq!(candidate.coating().cells()[&id].mass_kg, 0.01);
        assert_eq!(
            coating.reservoir,
            VaporReservoir::Prescribed {
                partial_pressure_pa: gas
                    .species_partial_pressure_pa(&fixture::species())
                    .unwrap(),
                temperature_k: gas.state().state().temperature,
            }
        );
    }
}

#[test]
fn bare_substrates_gain_real_coatings_and_heat_while_one_gas_cools_or_warms() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 350.0, 6000.0);
    let bindings = [fixture::binding(0.0), fixture::binding(0.0)];
    let chart = fixture::chart(&catalog, &bindings[0]);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let substrate = Substrate { datum: -100000.0 };
    let law = fixture::phase();
    let original = initial(
        &substrate,
        &law,
        &deposit_h,
        &bindings,
        vec![
            birth(0.0, cell(1, substrate, 2.0)),
            birth(0.0, cell(2, substrate, 0.0)),
        ],
        [0.0; 2],
        16,
    );
    let selection = selected(&bindings);
    let accepted = advance_coatings_with_gas(
        &original,
        &selection,
        None,
        &gas,
        &catalog,
        0.002,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let replay = advance_coatings_with_gas(
        &original,
        &[selection[1], selection[0]],
        None,
        &gas,
        &catalog,
        0.002,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(accepted.coating().cells(), replay.coating().cells());
    assert_eq!(accepted.coating().coatings(), replay.coating().coatings());
    assert_eq!(accepted.gas(), replay.gas());
    assert_eq!(accepted.audit(), replay.audit());
    for &id in &[CellId(1), CellId(2)] {
        assert!(accepted.coating().coatings()[&id].mass_kg > 0.0);
        assert!(accepted.coating().temperature_k(id).unwrap() > 300.0);
    }
    assert_ne!(accepted.gas().state().state().temperature, 350.0);
    close(
        accepted
            .gas()
            .species_mass_kg(&SpeciesId::new("A").unwrap())
            .unwrap(),
        gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap(),
        1e-14,
    );
    close(
        accepted.coating().deposit_mass_kg() + accepted.gas().mass_kg(),
        gas.mass_kg(),
        1e-14,
    );
    close(
        accepted.audit().solid_energy_change_j + accepted.audit().gas_energy_change_j
            - accepted.audit().incoming_energy_j
            - accepted.audit().thermal.external_heat_j,
        0.0,
        1e-7,
    );
    assert_eq!(accepted.audit().phase_steps, 2);
    assert!(original.cells().is_empty());
    check(&accepted);
}

#[test]
fn endpoint_birth_and_contact_have_no_retroactive_exposure_and_resume_exactly() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 350.0, 6000.0);
    let bindings = [fixture::binding(0.0), fixture::binding(0.0)];
    let chart = fixture::chart(&catalog, &bindings[0]);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let substrate = Substrate { datum: 0.0 };
    let law = fixture::phase();
    let events = vec![
        birth(0.0, cell(1, substrate, 0.0)),
        birth(0.001, cell(2, substrate, 0.0)),
        TimedEvent {
            time_s: 0.001,
            event: Event::Contact(ThermalLink {
                id: 9,
                a: CellId(1),
                b: CellId(2),
                conductance_w_k: 0.0,
                area_m2: 1e-4,
            }),
        },
    ];
    let original = initial(
        &substrate, &law, &deposit_h, &bindings, events, [0.0; 2], 16,
    );
    let selection = selected(&bindings);
    let first = advance_coatings_with_gas(
        &original,
        &selection,
        None,
        &gas,
        &catalog,
        0.001,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert!(first.coating().coatings()[&CellId(1)].mass_kg > 0.0);
    assert_eq!(first.coating().coatings()[&CellId(2)].mass_kg, 0.0);
    assert_eq!(first.audit().cells[&CellId(2)].phase_steps, 0);
    check(&first);
    let continued = advance_coatings_with_gas(
        first.coating(),
        &selection,
        None,
        first.gas(),
        &catalog,
        0.002,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let whole = advance_coatings_with_gas(
        &original,
        &selection,
        None,
        &gas,
        &catalog,
        0.002,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(continued.coating().cells(), whole.coating().cells());
    assert_eq!(continued.coating().coatings(), whole.coating().coatings());
    assert_eq!(continued.gas(), whole.gas());
    assert_eq!(continued.gas(), first.gas());
    assert_eq!(continued.coating().coatings()[&CellId(2)].mass_kg, 0.0);
    check(&continued);
    check(&whole);
}

#[test]
fn original_coating_heat_capacity_and_carrier_accommodation_share_the_same_storage() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 350.0, 6000.0);
    let mut bindings = [fixture::binding(0.0), fixture::binding(0.0)];
    for binding in &mut bindings {
        binding.accommodation_coefficient = 0.0;
    }
    let chart = fixture::chart(&catalog, &bindings[0]);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let substrate = Substrate { datum: 0.0 };
    let law = fixture::phase();
    let original = initial(
        &substrate,
        &law,
        &deposit_h,
        &bindings,
        vec![
            birth(0.0, cell(1, substrate, 4.0)),
            birth(0.0, cell(2, substrate, 4.0)),
        ],
        [0.002, 0.0],
        16,
    );
    let selection = selected(&bindings);
    let heated = advance_coatings_with_gas(
        &original,
        &selection,
        None,
        &gas,
        &catalog,
        0.01,
        heat(0.01),
        controls(),
        || Ok(()),
    )
    .unwrap();
    close(
        heated.coating().temperature_k(CellId(1)).unwrap(),
        300.0 + 0.04 / (8.0 + 0.002 * chart.1),
        1e-8,
    );
    close(
        heated.coating().temperature_k(CellId(2)).unwrap(),
        300.0 + 0.04 / 8.0,
        1e-8,
    );
    assert_eq!(heated.gas(), &gas);
    let accommodation = ThermalAccommodation::new(
        &catalog,
        &[
            (SpeciesId::new("A").unwrap(), 0.4),
            (fixture::species(), 0.3),
        ],
        [200.0, 900.0],
    )
    .unwrap();
    let maps = BTreeMap::from([
        (CellId(1), accommodation.clone()),
        (CellId(2), accommodation),
    ]);
    let reflected = advance_coatings_with_gas(
        &original,
        &selection,
        Some(&maps),
        &gas,
        &catalog,
        0.001,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(reflected.audit().transferred_mass_kg, 0.0);
    assert!(reflected.audit().reflected_heat_into_gas_j < 0.0);
    assert!(reflected.gas().state().state().temperature < gas.state().state().temperature);
    assert_eq!(reflected.coating().coatings()[&CellId(1)].mass_kg, 0.002);
    assert_eq!(reflected.coating().coatings()[&CellId(2)].mass_kg, 0.0);
    close(
        reflected.audit().energy_into_gas_j,
        reflected.audit().reflected_heat_into_gas_j,
        1e-10,
    );
    check(&heated);
    check(&reflected);
}

#[test]
fn independent_substrate_and_species_datums_preserve_small_physical_changes() {
    let base = fixture::catalog(0.0);
    let chart = fixture::chart(&base, &fixture::binding(0.0));
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let mut reference = None;
    for (datum, formation) in [(0.0, 0.0), (1e9, 1e12)] {
        let catalog = fixture::catalog(formation);
        let gas = fixture::gas(&catalog, 350.0, 6000.0);
        let bindings = [fixture::binding(formation), fixture::binding(formation)];
        let substrate = Substrate { datum };
        let original = initial(
            &substrate,
            &law,
            &deposit_h,
            &bindings,
            vec![
                birth(0.0, cell(1, substrate, 0.0)),
                birth(0.001, cell(2, substrate, 0.0)),
            ],
            [0.0, 1e-5],
            16,
        );
        let mut solve = controls();
        solve.chamber.surface.latent_compatibility_tolerance_j_kg = 0.01;
        solve.storage.absolute_energy_tolerance_j = 1e-8;
        let accepted = advance_coatings_with_gas(
            &original,
            &selected(&bindings),
            None,
            &gas,
            &catalog,
            0.002,
            heat(0.001),
            solve,
            || Ok(()),
        )
        .unwrap();
        let physical = (
            accepted.coating().temperature_k(CellId(1)).unwrap(),
            accepted.coating().temperature_k(CellId(2)).unwrap(),
            accepted.gas().state().state().temperature,
            accepted.coating().deposit_mass_kg(),
        );
        if let Some((a, b, t, mass)) = reference {
            close(physical.0, a, 1e-7);
            close(physical.1, b, 1e-7);
            close(physical.2, t, 1e-7);
            close(physical.3, mass, 1e-14);
        } else {
            reference = Some(physical);
        }
        check(&accepted);
    }
}

#[test]
fn cumulative_budget_cancellation_and_bad_kinetics_never_publish_partial_states() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 350.0, 6000.0);
    let bindings = [fixture::binding(0.0), fixture::binding(0.0)];
    let chart = fixture::chart(&catalog, &bindings[0]);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let substrate = Substrate { datum: 0.0 };
    let original = initial(
        &substrate,
        &law,
        &deposit_h,
        &bindings,
        vec![
            birth(0.0, cell(1, substrate, 0.0)),
            birth(0.0, cell(2, substrate, 0.0)),
        ],
        [0.0; 2],
        1,
    );
    let selection = selected(&bindings);
    let result = advance_coatings_with_gas(
        &original,
        &selection,
        None,
        &gas,
        &catalog,
        0.002,
        heat(0.001),
        controls(),
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(GasCoatingTimelineError::Thermal(
            ThermalSublimationError::Thermal(DepositionError::Budget(_))
        ))
    ));
    let mut polls = 0;
    let _accepted = advance_coatings_with_gas(
        &original,
        &selection,
        None,
        &gas,
        &catalog,
        0.001,
        heat(0.001),
        controls(),
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        let refused = advance_coatings_with_gas(
            &original,
            &selection,
            None,
            &gas,
            &catalog,
            0.001,
            heat(0.001),
            controls(),
            || {
                count += 1;
                if count == stop {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        );
        assert!(matches!(
            refused,
            Err(GasCoatingTimelineError::Thermal(
                ThermalSublimationError::Thermal(DepositionError::Cancelled)
            ))
        ));
    }
    let bad = fixture::binding(0.0);
    let mut wrong = bad.clone();
    wrong.accommodation_coefficient = 0.9;
    let refused = advance_coatings_with_gas(
        &original,
        &[
            GasCoatingTimelineSurface {
                cell_id: CellId(1),
                binding: &wrong,
            },
            selection[1],
        ],
        None,
        &gas,
        &catalog,
        0.001,
        heat(0.001),
        controls(),
        || Ok(()),
    );
    assert!(matches!(refused, Err(GasCoatingTimelineError::Thermal(_))));
    assert!(original.cells().is_empty());
    assert_eq!(original.pending_events().len(), 2);
    assert_eq!(original.time_s(), 0.0);
}

#[test]
fn zero_windows_and_endpoint_births_admit_actual_species_and_deposit_datums() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 350.0, 6000.0);
    let bindings = [fixture::binding(0.0), fixture::binding(0.0)];
    let chart = fixture::chart(&catalog, &bindings[0]);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let substrate = Substrate { datum: 0.0 };
    for time in [0.0, 0.5] {
        let original = initial(
            &substrate,
            &law,
            &deposit_h,
            &bindings,
            vec![
                birth(time, cell(1, substrate, 0.0)),
                birth(time, cell(2, substrate, 0.0)),
            ],
            [1e-5; 2],
            1,
        );
        let accepted = advance_coatings_with_gas(
            &original,
            &selected(&bindings),
            None,
            &gas,
            &catalog,
            time,
            heat(0.001),
            controls(),
            || Ok(()),
        )
        .unwrap();
        assert_eq!(accepted.audit().phase_steps, 0);
        assert_eq!(accepted.coating().deposit_mass_kg(), 2e-5);
        assert_eq!(accepted.gas(), &gas);
        let wrong = [fixture::binding(100.0), fixture::binding(0.0)];
        let refused = advance_coatings_with_gas(
            &original,
            &selected(&wrong),
            None,
            &gas,
            &catalog,
            time,
            heat(0.001),
            controls(),
            || Ok(()),
        );
        assert!(matches!(
            refused,
            Err(GasCoatingTimelineError::Exchange(
                GasSurfaceError::CaloricMismatch { .. }
            ))
        ));
        assert!(original.cells().is_empty());
        check(&accepted);
    }
}
