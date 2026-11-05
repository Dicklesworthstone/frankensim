use super::*;
use crate::enthalpy::{
    Boundary, CellGeometry, DepositionSimulation, Event, MaterialId, ThermalLaw, ThermalLink,
    TimedEvent,
};
use crate::sublimation::{
    SublimationError, SublimationLaw, SublimationProperties, SurfaceKinetics,
};

const SUBSTRATE_MASS: f64 = 0.02;
const COATING_MASS: f64 = 0.005;
const AREA: f64 = 0.001;

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = h / 500.0;
        if (100.0..=400.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 500.0
    }
    fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        // A source-specific admission stricter than its temperature chart.
        if cell.specific_h_j_kg > 375.0 * 500.0 {
            Err(DepositionError::MaterialDomain)
        } else {
            Ok(())
        }
    }
}
static CALORIC: Caloric = Caloric;

struct Phase;
impl SublimationLaw for Phase {
    fn temperature_range_k(&self) -> [f64; 2] {
        [100.0, 400.0]
    }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        let pressure = 10.0 * (t / 230.0).sqrt();
        Ok(SublimationProperties {
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: pressure / (2.0 * t),
            latent_heat_j_kg: 2.8e6,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 900.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}
static PHASE: Phase = Phase;
fn deposit_h(t: f64) -> Result<f64, SublimationError> {
    Ok(2000.0 * t)
}
fn bath() -> VaporReservoir {
    VaporReservoir::Prescribed {
        partial_pressure_pa: 20.0,
        temperature_k: 230.0,
    }
}
fn finite_bath() -> VaporReservoir {
    VaporReservoir::Finite {
        vapor_mass_kg: 0.001,
        volume_m3: 1.0,
        temperature_k: 230.0,
    }
}
fn cell(id: u64, t: f64, source_w: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: 7,
        material_id: MaterialId(1),
        mass_kg: SUBSTRATE_MASS,
        specific_h_j_kg: 500.0 * t,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 0.01, 0.0],
            volume_m3: 1e-5,
            surface_area_m2: AREA,
        },
        boundary: Boundary {
            ambient_k: 230.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w,
        },
    }
}
fn birth(time_s: f64, cell: DepositedCell) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Birth(cell),
    }
}
fn surface(id: u64, mass: f64, reservoir: VaporReservoir) -> CoatingSurface<'static> {
    CoatingSurface {
        cell_id: CellId(id),
        law: &PHASE,
        deposit_enthalpy: &deposit_h,
        kinetics: SurfaceKinetics {
            accommodation_coefficient: 1.0,
            mass_transfer_coefficient_m_s: None,
        },
        initial_mass_kg: mass,
        reservoir,
    }
}
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-6,
    }
}
fn solver() -> SolverConfig {
    SolverConfig {
        relative_tolerance: 1e-11,
        absolute_energy_tolerance_j: 1e-9,
        max_iterations: 128,
    }
}
fn thermal(events: Vec<TimedEvent>, steps: usize) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &CALORIC)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 12,
            max_substeps: steps,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0
}
fn simulation(
    events: Vec<TimedEvent>,
    surfaces: &[CoatingSurface<'static>],
    steps: usize,
) -> CoatingSimulation<'static> {
    CoatingSimulation::new(thermal(events, steps), surfaces, || Ok(())).unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-7, "{a:e} != {b:e}");
}
fn unchanged(left: &CoatingSimulation<'_>, right: &CoatingSimulation<'_>) {
    assert_eq!(left.time_s(), right.time_s());
    assert_eq!(left.cells(), right.cells());
    assert_eq!(left.coatings(), right.coatings());
    assert_eq!(left.links(), right.links());
    assert_eq!(left.pending_events(), right.pending_events());
}
fn exchange(
    state: &CoatingSimulation<'_>,
    id: u64,
    mass: f64,
    temperature: f64,
) -> CoatingCellExchange {
    let id = CellId(id);
    let expected = state.cells()[&id].clone();
    CoatingCellExchange {
        replacement: SupportedSurfaceState {
            substrate_mass_kg: expected.mass_kg,
            substrate_specific_enthalpy_j_kg: 500.0 * temperature,
            deposit_mass_kg: mass,
            area_m2: AREA,
        },
        expected_substrate: expected,
        expected_coating: state.coatings()[&id],
        interval_exposed_area_m2: AREA,
        replacement_reservoir: bath(),
    }
}

#[test]
fn external_heat_keeps_coating_capacity_and_never_runs_phantom_phase_exchange() {
    let link = ThermalLink {
        id: 9,
        a: CellId(1),
        b: CellId(2),
        area_m2: 0.0002,
        conductance_w_k: 0.5,
    };
    let state = simulation(
        vec![
            birth(0.0, cell(1, 230.0, 10.0)),
            birth(0.0, cell(2, 250.0, 0.0)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(link),
            },
        ],
        &[surface(1, COATING_MASS, bath()), surface(2, 0.0, bath())],
        8,
    );
    let (next, audit) = state
        .advance_heat_to_limited(0.25, config(0.25), solver(), 8, || Ok(()))
        .unwrap();
    // Cell 1 receives 10 W of conduction and 10 W of source power. Its
    // coating doubles its 10 J/K substrate capacity. Bare cell 2 loses 10 W.
    close(next.temperature_k(CellId(1)).unwrap(), 230.25);
    close(next.temperature_k(CellId(2)).unwrap(), 249.75);
    close(
        next.coatings()[&CellId(1)].specific_enthalpy_j_kg,
        2000.0 * 230.25,
    );
    assert_eq!(next.coatings()[&CellId(1)].mass_kg, COATING_MASS);
    assert_eq!(next.coatings()[&CellId(2)].mass_kg, 0.0);
    assert_eq!(next.mass_kg(), state.mass_kg());
    assert_eq!(next.coatings()[&CellId(1)].reservoir, bath());
    close(audit.stored_energy_change_j, 2.5);
    close(audit.external_heat_j, 2.5);
    close(audit.residual_j, 0.0);
    assert_eq!(audit.substeps, 1);
    assert_eq!(next.links()[&9], link);
    let exposed = next
        .exposed_areas_m2(&[CellId(2), CellId(1)], || Ok(()))
        .unwrap();
    assert_eq!(exposed[&CellId(1)], AREA - link.area_m2);
    assert_eq!(next.work_budget().max_substeps, 8);
    assert_eq!(next.coating_surfaces().len(), 2);
    assert_eq!(
        next.temperature_at_specific_enthalpy_k(CellId(1), 125000.0)
            .unwrap(),
        250.0
    );
}

#[test]
fn external_endpoints_grow_exhaust_and_regrow_only_the_deposit() {
    let original = cell(1, 230.0, 0.0);
    let mut state = simulation(
        vec![birth(0.0, original.clone())],
        &[surface(1, 0.0, bath())],
        8,
    );
    for (mass, temperature) in [(COATING_MASS, 240.0), (0.0, 235.0), (0.001, 245.0)] {
        let update = exchange(&state, 1, mass, temperature);
        state.apply_coating_exchanges(&[update], || Ok(())).unwrap();
        assert_eq!(state.coatings()[&CellId(1)].mass_kg, mass);
        assert_eq!(
            state.coatings()[&CellId(1)].specific_enthalpy_j_kg,
            2000.0 * temperature
        );
        assert_eq!(state.temperature_k(CellId(1)).unwrap(), temperature);
        let mut retained = state.cells()[&CellId(1)].clone();
        retained.specific_h_j_kg = original.specific_h_j_kg;
        assert_eq!(retained, original);
        assert_eq!(state.time_s(), 0.0);
    }
}

#[test]
fn heat_only_births_enter_at_the_endpoint_and_use_one_remaining_substep_allowance() {
    let state = simulation(
        vec![birth(0.5, cell(1, 230.0, 20.0))],
        &[surface(1, COATING_MASS, bath())],
        4,
    );
    let (born, audit) = state
        .advance_heat_to_limited(0.5, config(0.25), solver(), 0, || Ok(()))
        .unwrap();
    assert_eq!(audit.substeps, 0);
    assert_eq!(audit.events, 1);
    assert_eq!(audit.incoming_mass_kg, SUBSTRATE_MASS + COATING_MASS);
    close(
        audit.incoming_energy_j,
        SUBSTRATE_MASS * 500.0 * 230.0 + COATING_MASS * 2000.0 * 230.0,
    );
    assert_eq!(born.temperature_k(CellId(1)).unwrap(), 230.0);
    assert_eq!(born.coatings()[&CellId(1)].mass_kg, COATING_MASS);
    let (heated, heat) = born
        .advance_heat_to_limited(0.75, config(0.25), solver(), 1, || Ok(()))
        .unwrap();
    close(heated.temperature_k(CellId(1)).unwrap(), 230.25);
    assert_eq!(heat.substeps, 1);
    let before = heated.clone();
    assert!(matches!(
        heated.advance_heat_to_limited(1.0, config(0.25), solver(), 0, || Ok(())),
        Err(ThermalSublimationError::Thermal(DepositionError::Budget(_)))
    ));
    unchanged(&heated, &before);
    // A caller cannot enlarge the original four-substep ceiling either.
    assert!(matches!(
        born.advance_heat_to_limited(2.0, config(0.25), solver(), usize::MAX, || Ok(())),
        Err(ThermalSublimationError::Thermal(DepositionError::Budget(_)))
    ));
}

#[test]
fn stale_or_invalid_later_endpoints_refuse_the_complete_supported_batch() {
    let state = simulation(
        vec![
            birth(0.0, cell(1, 230.0, 0.0)),
            birth(0.0, cell(2, 230.0, 0.0)),
        ],
        &[surface(1, 0.0, bath()), surface(2, COATING_MASS, bath())],
        8,
    );
    let first = exchange(&state, 1, COATING_MASS, 240.0);
    let valid_second = exchange(&state, 2, 0.0, 235.0);
    let mut invalid = Vec::new();
    let mut stale_coating = valid_second.clone();
    stale_coating.expected_coating.mass_kg *= 0.5;
    invalid.push(stale_coating);
    let mut stale_substrate = valid_second.clone();
    stale_substrate.expected_substrate.boundary.source_w = 1.0;
    invalid.push(stale_substrate);
    let mut consumed_substrate = valid_second.clone();
    consumed_substrate.replacement.substrate_mass_kg *= 0.5;
    invalid.push(consumed_substrate);
    // This T is inside both charts but violates original validate_cell.
    invalid.push(exchange(&state, 2, 0.0, 390.0));
    let mut finite_stock = valid_second.clone();
    finite_stock.replacement_reservoir = finite_bath();
    invalid.push(finite_stock);
    for second in invalid {
        let mut trial = state.clone();
        assert!(
            trial
                .apply_coating_exchanges(&[first.clone(), second], || Ok(()))
                .is_err()
        );
        unchanged(&trial, &state);
    }
    let rows = [first, valid_second];
    let mut calls = 0;
    state
        .clone()
        .apply_coating_exchanges(&rows, || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut trial = state.clone();
        let mut current = 0;
        assert_eq!(
            trial.apply_coating_exchanges(&rows, || {
                current += 1;
                if current == stop {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(ThermalSublimationError::Thermal(DepositionError::Cancelled))
        );
        unchanged(&trial, &state);
    }
}

#[test]
fn external_coupling_refuses_private_future_and_shared_finite_gas_stocks() {
    for time in [0.0, 0.5] {
        let mut state = simulation(
            vec![birth(time, cell(1, 230.0, 0.0))],
            &[surface(1, 0.0, finite_bath())],
            8,
        );
        let before = state.clone();
        assert!(
            state
                .advance_heat_to_limited(1.0, config(1.0), solver(), 8, || Ok(()))
                .is_err()
        );
        assert!(state.apply_coating_exchanges(&[], || Ok(())).is_err());
        unchanged(&state, &before);
    }
    let mut shared = CoatingSimulation::new_shared_chamber(
        thermal(vec![birth(0.0, cell(1, 230.0, 0.0))], 8),
        &[super::super::ChamberCoatingSurface {
            cell_id: CellId(1),
            kinetics: surface(1, 0.0, bath()).kinetics,
            initial_mass_kg: 0.0,
        }],
        &PHASE,
        &deposit_h,
        finite_bath(),
        || Ok(()),
    )
    .unwrap();
    let before = shared.clone();
    assert!(
        shared
            .advance_heat_to_limited(1.0, config(1.0), solver(), 8, || Ok(()))
            .is_err()
    );
    assert!(shared.apply_coating_exchanges(&[], || Ok(())).is_err());
    unchanged(&shared, &before);
}
