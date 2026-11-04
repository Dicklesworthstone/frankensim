//! G0/G1: native heat storage, deposition inventories, events and rollback.

use super::{CoatingAudit, CoatingSimulation, CoatingSurface};
use crate::enthalpy::sublimation::ThermalSublimationError;
use crate::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, Reservoir, StepConfig, ThermalLaw, ThermalLink,
    TimedEvent,
};
use crate::sublimation::{
    MOLAR_GAS_CONSTANT, SolverConfig, SublimationError, SublimationLaw, SublimationProperties,
    SurfaceKinetics, VaporReservoir,
};

const SUBSTRATE_MASS: f64 = 0.02;
const SUBSTRATE_CP: f64 = 500.0;
const DEPOSIT_CP: f64 = 2000.0;
const AREA: f64 = 0.001;
const LATENT: f64 = 2.8e6;
const INITIAL_T: f64 = 230.0;

struct Caloric {
    datum: f64,
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = (h - self.datum) / SUBSTRATE_CP;
        if (100.0..=400.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / SUBSTRATE_CP
    }
}
static CALORIC: Caloric = Caloric { datum: 0.0 };

// p_sat / sqrt(T) is constant: prescribed-pressure molecular transfer has an
// independent analytical value even when latent heat changes the endpoint T.
// This manufactured fixture tests coupling, not a material correlation.
struct ManufacturedLaw;
impl SublimationLaw for ManufacturedLaw {
    fn temperature_range_k(&self) -> [f64; 2] {
        [100.0, 400.0]
    }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        let p = 10.0 * (t / INITIAL_T).sqrt();
        Ok(SublimationProperties {
            saturation_pressure_pa: p,
            saturation_pressure_derivative_pa_k: p / (2.0 * t),
            latent_heat_j_kg: LATENT,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 900.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}
static PHASE: ManufacturedLaw = ManufacturedLaw;

fn deposit_h(t: f64) -> Result<f64, SublimationError> {
    if (100.0..=400.0).contains(&t) {
        Ok(DEPOSIT_CP * t)
    } else {
        Err(SublimationError::Invalid("test deposit caloric domain"))
    }
}
fn bath(pressure_pa: f64) -> VaporReservoir {
    VaporReservoir::Prescribed {
        partial_pressure_pa: pressure_pa,
        temperature_k: INITIAL_T,
    }
}
fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: INITIAL_T,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}
fn cell(id: u64, t: f64, source_w: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id + 10,
        body_id: 7,
        material_id: MaterialId(1),
        mass_kg: SUBSTRATE_MASS,
        specific_h_j_kg: SUBSTRATE_CP * t,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 0.01, 0.0],
            volume_m3: 1e-5,
            surface_area_m2: AREA,
        },
        boundary: boundary(source_w),
    }
}
fn birth(time_s: f64, cell: DepositedCell) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Birth(cell),
    }
}
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-6,
    }
}
fn phase() -> SolverConfig {
    SolverConfig {
        relative_tolerance: 1e-11,
        absolute_energy_tolerance_j: 1e-9,
        max_iterations: 128,
    }
}
fn surface<'a>(
    id: u64,
    mass: f64,
    accommodation: f64,
    reservoir: VaporReservoir,
    enthalpy: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
) -> CoatingSurface<'a> {
    CoatingSurface {
        cell_id: CellId(id),
        law: &PHASE,
        deposit_enthalpy: enthalpy,
        kinetics: SurfaceKinetics {
            accommodation_coefficient: accommodation,
            mass_transfer_coefficient_m_s: None,
        },
        initial_mass_kg: mass,
        reservoir,
    }
}
fn simulation<'a>(
    events: Vec<TimedEvent>,
    caloric: &'a dyn ThermalLaw,
    surfaces: &[CoatingSurface<'a>],
    max_steps: usize,
) -> CoatingSimulation<'a> {
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(1), caloric)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 12,
            max_substeps: max_steps,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0;
    CoatingSimulation::new(thermal, surfaces, || Ok(())).unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} (tolerance {tolerance:e})"
    );
}
fn closes(audit: &CoatingAudit) {
    close(audit.mass_balance_residual_kg, 0.0, 1e-14);
    close(audit.energy_residual_j, 0.0, 1e-6);
    assert!(audit.absolute_energy_residual_j <= 1e-6);
}
fn growth_per_second() -> f64 {
    10.0 * AREA * (0.018 / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT * INITIAL_T)).sqrt()
}
fn same_state(left: &CoatingSimulation<'_>, right: &CoatingSimulation<'_>) {
    assert_eq!(left.time_s().to_bits(), right.time_s().to_bits());
    assert_eq!(left.cells(), right.cells());
    assert_eq!(left.links(), right.links());
    assert_eq!(left.pending_events(), right.pending_events());
    assert_eq!(left.coatings().len(), right.coatings().len());
    for (id, coat) in left.coatings() {
        let other = &right.coatings()[id];
        assert_eq!(coat.mass_kg.to_bits(), other.mass_kg.to_bits());
        assert_eq!(
            coat.specific_enthalpy_j_kg.to_bits(),
            other.specific_enthalpy_j_kg.to_bits()
        );
        assert_eq!(coat.reservoir, other.reservoir);
    }
}

#[test]
fn native_conduction_and_source_heat_both_substrate_and_deposited_storage() {
    let coat_mass = 0.005;
    let link = ThermalLink {
        id: 9,
        a: CellId(1),
        b: CellId(2),
        area_m2: 0.0002,
        conductance_w_k: 0.5,
    };
    let state = simulation(
        vec![
            birth(0.0, cell(1, INITIAL_T, 10.0)),
            birth(0.0, cell(2, 250.0, 0.0)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(link),
            },
        ],
        &CALORIC,
        &[surface(1, coat_mass, 0.0, bath(0.0), &deposit_h)],
        8,
    );
    let (next, audit) = state
        .advance_to(0.25, config(0.25), phase(), || Ok(()))
        .unwrap();
    // 10 W from the warmer neighbor plus 10 W authored source. The coating
    // doubles cell 1's capacity, while cell 2 loses the same conductive joules.
    close(next.temperature_k(CellId(1)).unwrap(), 230.25, 1e-9);
    close(next.temperature_k(CellId(2)).unwrap(), 249.75, 1e-9);
    close(
        next.coatings()[&CellId(1)].specific_enthalpy_j_kg,
        DEPOSIT_CP * 230.25,
        1e-6,
    );
    close(next.energy_j() - state.energy_j(), 2.5, 2e-8);
    close(audit.thermal.external_heat_j, 2.5, 1e-14);
    close(audit.stored_energy_change_j, 2.5, 2e-8);
    assert_eq!(next.cells()[&CellId(1)].mass_kg, SUBSTRATE_MASS);
    assert_eq!(next.coatings()[&CellId(1)].mass_kg, coat_mass);
    assert_eq!(next.links()[&9], link);
    assert_eq!(audit.transferred_mass_kg, 0.0);
    closes(&audit);
}

#[test]
fn zero_inventory_grows_exhausts_and_regrows_without_consuming_the_substrate() {
    let original = cell(1, INITIAL_T, 0.0);
    let state = simulation(
        vec![birth(0.0, original.clone())],
        &CALORIC,
        &[surface(1, 0.0, 1.0, bath(20.0), &deposit_h)],
        16,
    );
    let growth = growth_per_second();
    let (mut grown, growth_audit) = state
        .advance_to(1.0, config(1.0), phase(), || Ok(()))
        .unwrap();
    close(grown.coatings()[&CellId(1)].mass_kg, growth, 1e-17);
    close(
        grown.temperature_k(CellId(1)).unwrap(),
        INITIAL_T + growth * LATENT / 10.0,
        1e-8,
    );
    assert!(growth_audit.vapor_energy_j < 0.0);
    closes(&growth_audit);
    grown
        .set_reservoir(CellId(1), bath(0.0), || Ok(()))
        .unwrap();
    let (empty, loss_audit) = grown
        .advance_to(3.0, config(2.0), phase(), || Ok(()))
        .unwrap();
    assert_eq!(empty.coatings()[&CellId(1)].mass_kg, 0.0);
    close(loss_audit.transferred_mass_kg, growth, 1e-17);
    let expected =
        grown.temperature_k(CellId(1)).unwrap() - growth * LATENT / (10.0 + growth * DEPOSIT_CP);
    close(empty.temperature_k(CellId(1)).unwrap(), expected, 1e-8);
    let mut retained = empty.cells()[&CellId(1)].clone();
    retained.specific_h_j_kg = original.specific_h_j_kg;
    assert_eq!(retained, original);
    closes(&loss_audit);
    let (mut bare, bare_audit) = empty
        .advance_to(4.0, config(1.0), phase(), || Ok(()))
        .unwrap();
    assert_eq!(bare.coatings()[&CellId(1)].mass_kg, 0.0);
    assert_eq!(
        bare.temperature_k(CellId(1)),
        empty.temperature_k(CellId(1))
    );
    assert_eq!(bare_audit.transferred_mass_kg, 0.0);
    bare.set_reservoir(CellId(1), bath(20.0), || Ok(()))
        .unwrap();
    let (regrown, regrowth_audit) = bare
        .advance_to(4.5, config(0.5), phase(), || Ok(()))
        .unwrap();
    close(regrown.coatings()[&CellId(1)].mass_kg, 0.5 * growth, 1e-17);
    assert_eq!(regrown.cells()[&CellId(1)].mass_kg, SUBSTRATE_MASS);
    closes(&regrowth_audit);
}

#[test]
fn finite_vapor_supply_closes_actual_mass_and_sensible_plus_latent_energy() {
    let volume = 0.001;
    let equilibrium_mass = 10.0 * volume * 0.018 / (MOLAR_GAS_CONSTANT * INITIAL_T);
    let reservoir = VaporReservoir::Finite {
        vapor_mass_kg: 2.0 * equilibrium_mass,
        volume_m3: volume,
        temperature_k: INITIAL_T,
    };
    let state = simulation(
        vec![birth(0.0, cell(1, INITIAL_T, 0.0))],
        &CALORIC,
        &[surface(1, 0.0, 1.0, reservoir, &deposit_h)],
        4,
    );
    let (next, audit) = state
        .advance_to(1e4, config(1e4), phase(), || Ok(()))
        .unwrap();
    let coating = &next.coatings()[&CellId(1)];
    let VaporReservoir::Finite { vapor_mass_kg, .. } = coating.reservoir else {
        panic!("lost finite reservoir")
    };
    close(coating.mass_kg, equilibrium_mass, 1e-18);
    close(
        vapor_mass_kg + next.mass_kg(),
        2.0 * equilibrium_mass + state.mass_kg(),
        1e-17,
    );
    close(
        next.temperature_k(CellId(1)).unwrap(),
        INITIAL_T + equilibrium_mass * LATENT / 10.0,
        1e-8,
    );
    close(
        next.energy_j() - state.energy_j() + audit.vapor_energy_j,
        0.0,
        1e-7,
    );
    close(audit.transferred_mass_kg, -coating.mass_kg, 1e-18);
    closes(&audit);
}

#[test]
fn future_birth_and_endpoint_coverage_preserve_preceding_exposure_and_source_work() {
    let mut covered = boundary(6.0);
    covered.reservoir = Some(Reservoir {
        area_m2: AREA,
        conductance_w_k: 0.0,
        temperature_k: INITIAL_T,
    });
    let state = simulation(
        vec![
            birth(1.0, cell(1, INITIAL_T, 2.0)),
            TimedEvent {
                time_s: 2.0,
                event: Event::Boundary(CellId(1), covered),
            },
            TimedEvent {
                time_s: 3.0,
                event: Event::Boundary(CellId(1), boundary(0.0)),
            },
        ],
        &CALORIC,
        &[surface(1, 0.0, 1.0, bath(20.0), &deposit_h)],
        8,
    );
    let (waiting, idle) = state
        .advance_to(0.5, config(0.1), phase(), || Ok(()))
        .unwrap();
    assert!(waiting.cells().is_empty());
    assert_eq!(waiting.mass_kg(), 0.0);
    assert_eq!(idle.phase_steps, 0);
    assert_eq!(idle.thermal.substeps, 0);
    let (born, birth_audit) = waiting
        .advance_to(1.0, config(1.0), phase(), || Ok(()))
        .unwrap();
    assert_eq!(born.coatings()[&CellId(1)].mass_kg, 0.0);
    assert_eq!(birth_audit.phase_steps, 0);
    close(birth_audit.thermal.incoming_mass_kg, SUBSTRATE_MASS, 1e-18);
    let (covered_state, exposed_audit) = born
        .advance_to(2.0, config(1.0), phase(), || Ok(()))
        .unwrap();
    let growth = growth_per_second();
    close(covered_state.coatings()[&CellId(1)].mass_kg, growth, 1e-17);
    close(exposed_audit.thermal.external_heat_j, 2.0, 1e-14);
    assert_eq!(covered_state.cells()[&CellId(1)].boundary, covered);
    let (uncovered, covered_audit) = covered_state
        .advance_to(3.0, config(1.0), phase(), || Ok(()))
        .unwrap();
    assert_eq!(
        uncovered.coatings()[&CellId(1)].mass_kg,
        covered_state.coatings()[&CellId(1)].mass_kg
    );
    close(covered_audit.thermal.external_heat_j, 6.0, 1e-14);
    close(
        uncovered.temperature_k(CellId(1)).unwrap()
            - covered_state.temperature_k(CellId(1)).unwrap(),
        6.0 / (10.0 + growth * DEPOSIT_CP),
        1e-8,
    );
    let (resumed, resumed_audit) = uncovered
        .advance_to(3.5, config(0.5), phase(), || Ok(()))
        .unwrap();
    close(resumed.coatings()[&CellId(1)].mass_kg, 1.5 * growth, 1e-17);
    for audit in [&birth_audit, &exposed_audit, &covered_audit, &resumed_audit] {
        closes(audit);
    }
}

#[test]
fn cancellation_work_refusal_and_restart_keep_one_atomic_physical_history() {
    let make = |max_steps| {
        simulation(
            vec![birth(0.0, cell(1, INITIAL_T, 3.0))],
            &CALORIC,
            &[surface(1, 0.001, 1.0, bath(20.0), &deposit_h)],
            max_steps,
        )
    };
    let state = make(32);
    let original_cells = state.cells().clone();
    let initial_mass = state.mass_kg();
    let initial_energy = state.energy_j();
    let mut polls = 0;
    let error = state
        .advance_to(0.5, config(0.125), phase(), || {
            polls += 1;
            if polls >= 12 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        })
        .err()
        .expect("cancelled candidate must refuse");
    assert_eq!(
        error,
        ThermalSublimationError::Thermal(DepositionError::Cancelled)
    );
    assert_eq!(state.cells(), &original_cells);
    assert_eq!(state.mass_kg(), initial_mass);
    assert_eq!(state.energy_j(), initial_energy);
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.coatings()[&CellId(1)].mass_kg, 0.001);
    assert_eq!(state.coatings()[&CellId(1)].reservoir, bath(20.0));
    let limited = make(2);
    let error = limited
        .advance_to(0.375, config(0.125), phase(), || Ok(()))
        .err()
        .expect("third thermal step must exceed the shared budget");
    assert!(matches!(
        error,
        ThermalSublimationError::Thermal(DepositionError::Budget(_))
    ));
    same_state(&state, &limited);
    assert!(
        state
            .advance_to(
                0.125,
                config(0.125),
                SolverConfig {
                    max_iterations: 1,
                    ..phase()
                },
                || Ok(())
            )
            .is_err()
    );
    let (whole, whole_audit) = state
        .advance_to(1.0, config(0.125), phase(), || Ok(()))
        .unwrap();
    let (half, first_audit) = state
        .advance_to(0.5, config(0.125), phase(), || Ok(()))
        .unwrap();
    let (resumed, second_audit) = half
        .advance_to(1.0, config(0.125), phase(), || Ok(()))
        .unwrap();
    same_state(&whole, &resumed);
    close(
        whole_audit.transferred_mass_kg,
        first_audit.transferred_mass_kg + second_audit.transferred_mass_kg,
        1e-17,
    );
    close(
        whole_audit.vapor_energy_j,
        first_audit.vapor_energy_j + second_audit.vapor_energy_j,
        1e-8,
    );
    closes(&whole_audit);
    closes(&first_audit);
    closes(&second_audit);
}

#[test]
fn independent_large_caloric_datums_preserve_heated_deposition_and_its_audit() {
    let substrate_shift = 1e8;
    let deposit_shift = -3e8;
    let shifted_caloric = Caloric {
        datum: substrate_shift,
    };
    let shifted_deposit = |t| deposit_h(t).map(|h| h + deposit_shift);
    let original_cell = cell(1, INITIAL_T, 4.0);
    let mut shifted_cell = original_cell.clone();
    shifted_cell.specific_h_j_kg += substrate_shift;
    let original = simulation(
        vec![birth(0.0, original_cell)],
        &CALORIC,
        &[surface(1, 0.00125, 1.0, bath(20.0), &deposit_h)],
        4,
    );
    let shifted = simulation(
        vec![birth(0.0, shifted_cell)],
        &shifted_caloric,
        &[surface(1, 0.00125, 1.0, bath(20.0), &shifted_deposit)],
        4,
    );
    let (original, audit) = original
        .advance_to(0.25, config(0.25), phase(), || Ok(()))
        .unwrap();
    let (shifted, shifted_audit) = shifted
        .advance_to(0.25, config(0.25), phase(), || Ok(()))
        .unwrap();
    close(
        shifted.temperature_k(CellId(1)).unwrap(),
        original.temperature_k(CellId(1)).unwrap(),
        1e-8,
    );
    close(
        shifted.coatings()[&CellId(1)].mass_kg,
        original.coatings()[&CellId(1)].mass_kg,
        1e-17,
    );
    close(
        shifted_audit.vapor_energy_j - audit.vapor_energy_j,
        audit.transferred_mass_kg * deposit_shift,
        1e-6,
    );
    close(
        shifted_audit.thermal.external_heat_j,
        audit.thermal.external_heat_j,
        1e-14,
    );
    closes(&audit);
    closes(&shifted_audit);
}

#[test]
fn preexisting_coating_on_a_future_birth_carries_its_own_mass_and_sensible_energy() {
    let coat_mass = 0.005;
    let birth_t = 240.0;
    let state = simulation(
        vec![birth(2.0, cell(1, birth_t, 0.0))],
        &CALORIC,
        &[surface(1, coat_mass, 0.0, bath(0.0), &deposit_h)],
        0,
    );
    assert_eq!(state.mass_kg(), 0.0);
    assert_eq!(state.energy_j(), 0.0);
    let (born, audit) = state
        .advance_to(2.0, config(0.1), phase(), || Ok(()))
        .unwrap();
    let incoming_energy = (SUBSTRATE_MASS * SUBSTRATE_CP + coat_mass * DEPOSIT_CP) * birth_t;
    close(born.mass_kg(), SUBSTRATE_MASS + coat_mass, 1e-18);
    close(born.energy_j(), incoming_energy, 1e-10);
    close(
        audit.thermal.incoming_mass_kg,
        SUBSTRATE_MASS + coat_mass,
        1e-18,
    );
    close(audit.thermal.incoming_energy_j, incoming_energy, 1e-10);
    close(audit.stored_energy_change_j, incoming_energy, 1e-10);
    assert_eq!(born.cells()[&CellId(1)].mass_kg, SUBSTRATE_MASS);
    assert_eq!(born.coatings()[&CellId(1)].mass_kg, coat_mass);
    assert_eq!(audit.thermal.substeps, 0);
    assert_eq!(audit.phase_steps, 0);
    assert_eq!(audit.thermal.events, 1);
    assert!(born.pending_events().is_empty());
    assert_eq!(state.pending_events().len(), 1);
    closes(&audit);
}
