//! Actual thermal timeline plus original molecular/gas owners; synthetic calorics.
use super::*;
use crate::sublimation::gas::tests as fixture;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositionBudget, MaterialId, Reservoir, ThermalLaw, ThermalLink,
    TimedEvent,
};
use fs_thermochem::SpeciesId;

#[derive(Clone, Copy)]
struct Caloric(f64, f64);
impl Caloric {
    fn h(self, temperature: f64) -> f64 { self.0 + self.1 * (temperature - 400.0) }
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 400.0 + (h - self.0) / self.1;
        if (200.0..=900.0).contains(&t) { Ok(t) } else { Err(DepositionError::MaterialDomain) }
    }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 / self.1 }
}
fn chart(catalog: &ThermochemicalGasModel, binding: &GasSurfaceBinding) -> Caloric {
    let (h, cp) = fixture::chart(catalog, binding);
    Caloric(h, cp)
}
fn controls() -> GasChamberSolve {
    GasChamberSolve { surface: fixture::controls(), maximum_surfaces: 2 }
}
fn heat(dt: f64) -> StepConfig { StepConfig { max_dt_s: dt, energy_tolerance_j: 1e-7 } }
fn boundary(power: f64, covered_area: f64) -> Boundary {
    Boundary {
        ambient_k: 350.0, htc_w_m2_k: 0.0, source_w: power,
        reservoir: Some(Reservoir {
            area_m2: covered_area, conductance_w_k: 0.0, temperature_k: 350.0,
        }),
    }
}
fn cell(id: u64, temperature: f64, caloric: Caloric, power: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id), segment_id: id, body_id: id, material_id: MaterialId(1),
        mass_kg: 0.01, specific_h_j_kg: caloric.h(temperature),
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0], end_m: [id as f64, 0.01, 0.0],
            volume_m3: 0.01 / 917.0, surface_area_m2: 2e-4,
        },
        boundary: boundary(power, 1e-4),
    }
}
fn birth(time_s: f64, cell: DepositedCell) -> TimedEvent {
    TimedEvent { time_s, event: Event::Birth(cell) }
}
fn schedule(caloric: &Caloric, events: Vec<TimedEvent>, steps: usize) -> DepositionSimulation<'_> {
    DepositionSimulation::new(events, &[(MaterialId(1), caloric)], DepositionBudget {
        max_cells: 2, max_links: 1, max_events: 16, max_substeps: steps,
    }).unwrap()
}
fn selection<'a>(phase: &'a dyn SublimationLaw, binding: &'a GasSurfaceBinding,
    caloric: Caloric) -> [GasTimelineSurface<'a>; 2] {
    std::array::from_fn(|i| GasTimelineSurface {
        cell_id: CellId(i as u64 + 1), law: phase, binding,
        pure_solid_enthalpy_range_j_kg: [caloric.h(200.0), caloric.h(900.0)],
    })
}
fn law(catalog: &ThermochemicalGasModel, alpha: f64) -> ThermalAccommodation {
    ThermalAccommodation::new(catalog, &[
        (SpeciesId::new("A").unwrap(), alpha), (fixture::species(), alpha),
    ], [200.0, 900.0]).unwrap()
}
fn laws(catalog: &ThermochemicalGasModel, alpha: f64) -> BTreeMap<CellId, ThermalAccommodation> {
    BTreeMap::from([(CellId(1), law(catalog, alpha)), (CellId(2), law(catalog, alpha))])
}

#[test]
fn actual_heat_contact_and_phase_share_reflected_energy_without_a_second_heat_source() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let caloric = chart(&catalog, &binding);
    let phase = fixture::phase();
    let selected = selection(&phase, &binding, caloric);
    let thermal = schedule(&caloric, vec![
        birth(0.0, cell(1, 420.0, caloric, 2.0)),
        birth(0.0, cell(2, 300.0, caloric, 2.0)),
        TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id: 9, a: CellId(1), b: CellId(2), area_m2: 1e-5, conductance_w_k: 0.02,
        }) },
    ], 1000).advance_to(0.0, heat(1.0), || Ok(())).unwrap().0;
    let before = thermal.cells().clone();
    let gas = fixture::gas(&catalog, 350.0, 800.0);
    let coefficients = laws(&catalog, 0.5);
    let accepted = advance_deposition_with_accommodation(&thermal, &selected, &coefficients,
        &gas, &catalog, 0.03125, heat(0.00390625), controls(), || Ok(())).unwrap();
    let report = accepted.audit();
    assert!(report.cells[&CellId(1)].reflected_heat_into_gas_j > 0.0);
    assert!(report.cells[&CellId(2)].reflected_heat_into_gas_j < 0.0);
    assert!(report.cells[&CellId(1)].transferred_mass_kg > 0.0);
    assert!(report.cells[&CellId(2)].transferred_mass_kg < 0.0);
    fixture::near(report.thermal.external_heat_j, 0.125, 1e-12);
    fixture::near(report.reflected_heat_into_gas_j,
        report.cells.values().map(|row| row.reflected_heat_into_gas_j).sum(), 1e-10);
    fixture::near(report.solid_energy_change_j + report.gas_energy_change_j
        - report.incoming_energy_j - report.thermal.external_heat_j, 0.0, 1e-7);
    fixture::near(report.energy_into_gas_j, report.gas_energy_change_j, 1e-7);
    fixture::near(report.mass_residual_kg, 0.0, 1e-14);
    fixture::near(accepted.gas().species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap(),
        gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap(), 1e-14);
    assert_eq!(accepted.gas().point().axes().get("local-lot"), Some(&9.0));
    assert_eq!(accepted.thermal().links(), thermal.links());
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn birth_and_endpoint_cover_apply_to_reflection_at_the_same_physical_times() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let caloric = chart(&catalog, &binding);
    let phase = fixture::phase();
    let selected = selection(&phase, &binding, caloric);
    let thermal = schedule(&caloric, vec![
        birth(0.03125, cell(1, 400.0, caloric, 0.0)),
        birth(0.03125, cell(2, 340.0, caloric, 0.0)),
        TimedEvent { time_s: 0.0625, event: Event::Boundary(CellId(1), boundary(0.0, 2e-4)) },
        TimedEvent { time_s: 0.0625, event: Event::Boundary(CellId(2), boundary(0.0, 2e-4)) },
    ], 100);
    let gas = fixture::gas(&catalog, 300.0, 0.0);
    let coefficients = laws(&catalog, 0.5);
    let at_birth = advance_deposition_with_accommodation(&thermal, &selected, &coefficients,
        &gas, &catalog, 0.03125, heat(0.03125), controls(), || Ok(())).unwrap();
    assert_eq!(at_birth.gas(), &gas);
    assert_eq!(at_birth.audit().reflected_heat_into_gas_j, 0.0);
    assert_eq!(at_birth.audit().phase_steps, 0);
    let complete = advance_deposition_with_accommodation(&thermal, &selected, &coefficients,
        &gas, &catalog, 0.125, heat(0.03125), controls(), || Ok(())).unwrap();

    // Exactly one exposed interval exists. Compare directly with the same
    // chamber owner; no pre-birth or post-cover encounter is permitted.
    let invert = |h| fixture::temperature(h, (caloric.0, caloric.1));
    let rows = [400.0, 340.0].map(|t| fixture::solid(t, (caloric.0, caloric.1)));
    let surfaces: [GasChamberSurface<'_>; 2] = std::array::from_fn(|i| GasChamberSurface {
        id: i as u64 + 1, law: &phase, binding: &binding, solid: rows[i],
        temperature_from_enthalpy: &invert,
    });
    let chamber_laws = BTreeMap::from([(1, law(&catalog, 0.5)), (2, law(&catalog, 0.5))]);
    let (expected, reflected) = advance_chamber_with_accommodation(&surfaces, &chamber_laws,
        &gas, &catalog, 0.03125, controls(), || false).unwrap();
    assert_eq!(complete.gas(), expected.gas());
    for endpoint in expected.surfaces() {
        let id = CellId(endpoint.id);
        assert_eq!(complete.thermal().cells()[&id].mass_kg, endpoint.solid.mass_kg);
        assert_eq!(complete.thermal().cells()[&id].specific_h_j_kg, endpoint.solid.specific_enthalpy_j_kg);
        assert_eq!(complete.thermal().exposed_area_m2(id), Ok(0.0));
        fixture::near(complete.audit().cells[&id].reflected_heat_into_gas_j,
            reflected[&endpoint.id].heat_into_gas_j, 1e-10);
    }
    fixture::near(complete.audit().energy_residual_j, 0.0, 1e-7);
    assert!(thermal.cells().is_empty());
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn zero_accommodation_keeps_existing_timeline_bits_and_order_is_canonical() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let caloric = chart(&catalog, &binding);
    let phase = fixture::phase();
    let selected = selection(&phase, &binding, caloric);
    let thermal = schedule(&caloric, vec![birth(0.0, cell(1, 400.0, caloric, 1.0)),
        birth(0.0, cell(2, 300.0, caloric, 1.0))], 1000);
    let gas = fixture::gas(&catalog, 350.0, 800.0);
    let original = advance_deposition_with_gas(&thermal, &selected, &gas, &catalog,
        0.03125, heat(0.0078125), controls(), || Ok(())).unwrap();
    let zeros = advance_deposition_with_accommodation(&thermal, &selected, &laws(&catalog, 0.0),
        &gas, &catalog, 0.03125, heat(0.0078125), controls(), || Ok(())).unwrap();
    assert_eq!(zeros.gas(), original.gas());
    assert_eq!(zeros.thermal().cells(), original.thermal().cells());
    assert_eq!(zeros.audit(), original.audit());
    let coefficients = laws(&catalog, 0.3);
    let forward = advance_deposition_with_accommodation(&thermal, &selected, &coefficients,
        &gas, &catalog, 0.03125, heat(0.0078125), controls(), || Ok(())).unwrap();
    let reverse = advance_deposition_with_accommodation(&thermal, &[selected[1], selected[0]],
        &coefficients, &gas, &catalog, 0.03125, heat(0.0078125), controls(), || Ok(())).unwrap();
    assert_eq!(forward.gas(), reverse.gas());
    assert_eq!(forward.thermal().cells(), reverse.thermal().cells());
    assert_eq!(forward.audit(), reverse.audit());
    assert_ne!(forward.gas().state().state().temperature, original.gas().state().state().temperature);
}

#[test]
fn continuation_reuses_the_actual_endpoint_and_reflection_totals_are_additive() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let caloric = chart(&catalog, &binding);
    let phase = fixture::phase();
    let selected = selection(&phase, &binding, caloric);
    let thermal = schedule(&caloric, vec![birth(0.0, cell(1, 400.0, caloric, 0.5)),
        birth(0.0, cell(2, 300.0, caloric, 0.5))], 1000);
    let gas = fixture::gas(&catalog, 350.0, 800.0);
    let coefficients = laws(&catalog, 0.4);
    let first = advance_deposition_with_accommodation(&thermal, &selected, &coefficients,
        &gas, &catalog, 0.03125, heat(0.015625), controls(), || Ok(())).unwrap();
    let resumed = advance_deposition_with_accommodation(first.thermal(), &selected, &coefficients,
        first.gas(), &catalog, 0.0625, heat(0.015625), controls(), || Ok(())).unwrap();
    let whole = advance_deposition_with_accommodation(&thermal, &selected, &coefficients,
        &gas, &catalog, 0.0625, heat(0.015625), controls(), || Ok(())).unwrap();
    assert_eq!(resumed.gas(), whole.gas());
    assert_eq!(resumed.thermal().cells(), whole.thermal().cells());
    assert_eq!(resumed.audit().thermal.incoming_mass_kg, 0.0);
    fixture::near(first.audit().reflected_heat_into_gas_j + resumed.audit().reflected_heat_into_gas_j,
        whole.audit().reflected_heat_into_gas_j, 1e-10);
    for id in [CellId(1), CellId(2)] {
        fixture::near(first.audit().cells[&id].reflected_heat_into_gas_j
            + resumed.audit().cells[&id].reflected_heat_into_gas_j,
            whole.audit().cells[&id].reflected_heat_into_gas_j, 1e-10);
    }
}

#[test]
fn inactive_selected_cells_require_coefficients_but_cause_no_encounters_or_timestep_loop() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let caloric = chart(&catalog, &binding);
    let phase = fixture::phase();
    let selected = selection(&phase, &binding, caloric);
    let thermal = schedule(&caloric, vec![birth(1e9, cell(1, 400.0, caloric, 0.0)),
        birth(1e9, cell(2, 300.0, caloric, 0.0))], 1);
    let gas = fixture::gas(&catalog, 350.0, 800.0);
    let coefficients = laws(&catalog, 0.5);
    let mut calls = 0;
    let quiet = advance_deposition_with_accommodation(&thermal, &selected, &coefficients,
        &gas, &catalog, 1e8, heat(1e-8), controls(), || { calls += 1; Ok(()) }).unwrap();
    assert!(calls < 100);
    assert!(quiet.thermal().cells().is_empty());
    assert_eq!(quiet.gas(), &gas);
    assert_eq!(quiet.audit().reflected_heat_into_gas_j, 0.0);
    assert_eq!(quiet.audit().phase_steps, 0);
    for bad in [BTreeMap::new(), BTreeMap::from([(CellId(1), law(&catalog, 0.5))]),
        BTreeMap::from([(CellId(1), law(&catalog, 0.5)), (CellId(3), law(&catalog, 0.5))])] {
        assert!(advance_deposition_with_accommodation(&thermal, &selected, &bad,
            &gas, &catalog, 0.0, heat(0.01), controls(), || Ok(())).is_err());
    }
}

#[test]
fn reflected_caloric_exit_budget_and_late_cancellation_are_transactional() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let caloric = chart(&catalog, &binding);
    let phase = fixture::phase();
    let selected = selection(&phase, &binding, caloric);
    let events = vec![birth(0.0, cell(1, 400.0, caloric, 0.0)),
        birth(0.0, cell(2, 400.0, caloric, 0.0))];
    let thermal = schedule(&caloric, events.clone(), 1000)
        .advance_to(0.0, heat(1.0), || Ok(())).unwrap().0;
    let before = thermal.cells().clone();
    let gas = fixture::gas(&catalog, 300.0, 0.0);
    let saved = gas.clone();
    let coefficients = laws(&catalog, 0.5);
    let mut narrow = selected;
    narrow[1].pure_solid_enthalpy_range_j_kg[0] = caloric.h(399.999);
    assert!(matches!(advance_deposition_with_accommodation(&thermal, &narrow, &coefficients,
        &gas, &catalog, 0.01, heat(0.01), controls(), || Ok(())),
        Err(GasTimelineError::Thermal(DepositionError::MaterialDomain))));
    let limited = schedule(&caloric, events, 1);
    assert!(matches!(advance_deposition_with_accommodation(&limited, &selected, &coefficients,
        &gas, &catalog, 0.03125, heat(0.0078125), controls(), || Ok(())),
        Err(GasTimelineError::Thermal(DepositionError::Budget(_)))));
    let mut calls = 0;
    advance_deposition_with_accommodation(&thermal, &selected, &coefficients,
        &gas, &catalog, 0.03125, heat(0.0078125), controls(), || { calls += 1; Ok(()) }).unwrap();
    let mut count = 0;
    assert!(matches!(advance_deposition_with_accommodation(&thermal, &selected, &coefficients,
        &gas, &catalog, 0.03125, heat(0.0078125), controls(), || {
            count += 1;
            if count == calls { Err(DepositionError::Cancelled) } else { Ok(()) }
        }), Err(GasTimelineError::Thermal(DepositionError::Cancelled))));
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(gas, saved);
    assert!(limited.cells().is_empty());
}
