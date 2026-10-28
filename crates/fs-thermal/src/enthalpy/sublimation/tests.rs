//! G0 conservation/transactionality and G1 coupling-refinement checks.
use super::*;
use crate::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, MaterialId, Reservoir, ThermalLaw,
    ThermalLink, TimedEvent,
};
use crate::sublimation::SublimationProperties;

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 250.0 + h / 2000.0;
        if (180.0..=273.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 2000.0
    }
}
static CALORIC: Caloric = Caloric;

// Analytic Clausius-Clapeyron fixture: L = R_specific * 6000 K.
// It tests coupling, not the accuracy of any empirical water-ice correlation.
struct Ice;
impl SublimationLaw for Ice {
    fn temperature_range_k(&self) -> [f64; 2] {
        [180.0, 273.0]
    }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        if !(180.0..=273.0).contains(&t) {
            return Err(SublimationError::TemperatureOutOfRange {
                temperature_k: t,
                minimum_k: 180.0,
                maximum_k: 273.0,
            });
        }
        let p = 100.0 * (6000.0 * (1.0 / 250.0 - 1.0 / t)).exp();
        Ok(SublimationProperties {
            saturation_pressure_pa: p,
            saturation_pressure_derivative_pa_k: p * 6000.0 / (t * t),
            latent_heat_j_kg: 8.314_462_618_153_24 / 0.018 * 6000.0,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 917.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}
fn kinetics() -> SurfaceKinetics {
    SurfaceKinetics {
        accommodation_coefficient: 0.01,
        mass_transfer_coefficient_m_s: None,
    }
}
fn phase() -> SolverConfig {
    SolverConfig {
        relative_tolerance: 1e-11,
        absolute_energy_tolerance_j: 1e-9,
        max_iterations: 96,
    }
}
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-7,
    }
}
fn vacuum() -> VaporReservoir {
    VaporReservoir::Prescribed {
        partial_pressure_pa: 0.0,
        temperature_k: 250.0,
    }
}
fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 250.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}
fn cell(id: u64, h: f64, source_w: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: 1,
        material_id: MaterialId(1),
        mass_kg: 0.02,
        specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 0.01, 0.0],
            volume_m3: 0.02 / 917.0,
            surface_area_m2: 0.001,
        },
        boundary: boundary(source_w),
    }
}
fn scheduled(events: Vec<TimedEvent>, max_steps: usize) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &CALORIC)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: max_steps,
        },
    )
    .unwrap()
}
fn birth(time_s: f64, cell: DepositedCell) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Birth(cell),
    }
}
fn initial(source_w: f64) -> DepositionSimulation<'static> {
    scheduled(vec![birth(0.0, cell(1, 0.0, source_w))], 10_000)
        .advance_to(0.0, config(1.0), || Ok(()))
        .unwrap()
        .0
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn finite_vapor_and_actual_solid_share_mass_and_complete_energy_transfer() {
    let mut thermal = initial(0.0);
    let before = thermal.cells()[&CellId(1)].clone();
    let old_mass = thermal.mass_kg();
    let old_energy = thermal.energy_j();
    let mut reservoir = VaporReservoir::Finite {
        vapor_mass_kg: 0.0,
        volume_m3: 0.001,
        temperature_k: 250.0,
    };
    let result = thermal
        .sublimate_cell(
            CellId(1),
            &Ice,
            kinetics(),
            &mut reservoir,
            2.0,
            phase(),
            || Ok(()),
        )
        .unwrap();
    let VaporReservoir::Finite { vapor_mass_kg, .. } = reservoir else {
        panic!("finite vapor");
    };
    assert!(result.transferred_mass_kg > 0.0);
    assert!(thermal.temperature_k(CellId(1)).unwrap() < 250.0);
    close(thermal.mass_kg() + vapor_mass_kg, old_mass, 1e-17);
    close(thermal.energy_j() + result.vapor_energy_j, old_energy, 1e-7);
    close(
        result.vapor_energy_j - result.latent_energy_j,
        result.transferred_mass_kg * result.solid.specific_enthalpy_j_kg,
        1e-12,
    );
    assert_eq!(thermal.cells()[&CellId(1)].geometry, before.geometry);
    assert_eq!(thermal.time_s(), 0.0);
    assert!(result.recession_m > 0.0);
}

#[test]
fn deposition_grows_and_heats_the_actual_cell() {
    let mut thermal = initial(0.0);
    let old_mass = thermal.mass_kg();
    let mut reservoir = VaporReservoir::Prescribed {
        partial_pressure_pa: 200.0,
        temperature_k: 250.0,
    };
    let result = thermal
        .sublimate_cell(
            CellId(1),
            &Ice,
            kinetics(),
            &mut reservoir,
            1.0,
            phase(),
            || Ok(()),
        )
        .unwrap();
    assert!(result.transferred_mass_kg < 0.0);
    assert!(result.latent_energy_j < 0.0);
    assert!(result.recession_m < 0.0);
    assert!(thermal.mass_kg() > old_mass);
    assert!(thermal.temperature_k(CellId(1)).unwrap() > 250.0);
    close(thermal.energy_j() + result.vapor_energy_j, 0.0, 1e-7);
}

#[test]
fn covered_surface_exchanges_no_mass_or_latent_heat() {
    let mut covered = cell(1, 0.0, 0.0);
    covered.boundary.reservoir = Some(Reservoir {
        area_m2: covered.geometry.surface_area_m2,
        conductance_w_k: 0.0,
        temperature_k: 250.0,
    });
    let mut thermal = scheduled(vec![birth(0.0, covered)], 100)
        .advance_to(0.0, config(1.0), || Ok(()))
        .unwrap()
        .0;
    let before = thermal.cells().clone();
    let mut gas = vacuum();
    let result = thermal
        .sublimate_cell(CellId(1), &Ice, kinetics(), &mut gas, 1.0, phase(), || {
            Ok(())
        })
        .unwrap();
    assert_eq!(result.transferred_mass_kg, 0.0);
    assert_eq!(result.latent_energy_j, 0.0);
    assert_eq!(thermal.cells(), &before);
}

#[test]
fn conduction_source_and_sublimation_close_one_global_energy_account() {
    let thermal = scheduled(
        vec![
            birth(0.0, cell(1, 0.0, 2.0)),
            birth(0.0, cell(2, 20_000.0, 0.0)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 7,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.0001,
                    conductance_w_k: 0.4,
                }),
            },
        ],
        10_000,
    )
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0;
    let (next, _, audit) = thermal
        .advance_sublimating_to(
            CellId(1),
            &Ice,
            kinetics(),
            vacuum(),
            2.0,
            config(0.1),
            phase(),
            || Ok(()),
        )
        .unwrap();
    close(audit.thermal.external_heat_j, 4.0, 1e-12);
    close(
        next.energy_j() - thermal.energy_j() + audit.vapor_energy_j,
        4.0,
        1e-7,
    );
    close(
        next.mass_kg() - thermal.mass_kg() + audit.transferred_mass_kg,
        0.0,
        1e-16,
    );
    assert!(audit.thermal.substeps >= audit.phase_steps);
    assert!(next.temperature_k(CellId(2)).unwrap() < 260.0);
    assert_eq!(thermal.time_s(), 0.0);
    close(audit.energy_residual_j, 0.0, 1e-7);
}

#[test]
fn phase_refusal_and_cancellation_leave_both_owners_unchanged() {
    let mut thermal = initial(0.0);
    let original = thermal.cells().clone();
    let mut gas = VaporReservoir::Finite {
        vapor_mass_kg: 0.0,
        volume_m3: 0.001,
        temperature_k: 250.0,
    };
    let old_gas = gas.clone();
    let mut bad = phase();
    bad.max_iterations = 0;
    assert!(
        thermal
            .sublimate_cell(CellId(1), &Ice, kinetics(), &mut gas, 1.0, bad, || Ok(()))
            .is_err()
    );
    assert_eq!(thermal.cells(), &original);
    assert_eq!(gas, old_gas);
    let result =
        thermal.sublimate_cell(CellId(1), &Ice, kinetics(), &mut gas, 1.0, phase(), || {
            Err(DepositionError::Cancelled)
        });
    assert_eq!(
        result,
        Err(ThermalSublimationError::Thermal(DepositionError::Cancelled))
    );
    assert_eq!(thermal.cells(), &original);
    assert_eq!(gas, old_gas);
}

#[test]
fn endpoint_cover_event_does_not_remove_the_preceding_exposed_interval() {
    let mut covered = boundary(0.0);
    covered.reservoir = Some(Reservoir {
        area_m2: 0.001,
        conductance_w_k: 0.0,
        temperature_k: 250.0,
    });
    let thermal = scheduled(
        vec![
            birth(0.0, cell(1, 0.0, 0.0)),
            TimedEvent {
                time_s: 0.5,
                event: Event::Boundary(CellId(1), covered),
            },
        ],
        100,
    )
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0;
    let mut expected = initial(0.0);
    expected
        .sublimate_cell(
            CellId(1),
            &Ice,
            kinetics(),
            &mut vacuum(),
            0.5,
            phase(),
            || Ok(()),
        )
        .unwrap();
    let (next, _, audit) = thermal
        .advance_sublimating_to(
            CellId(1),
            &Ice,
            kinetics(),
            vacuum(),
            1.0,
            config(1.0),
            phase(),
            || Ok(()),
        )
        .unwrap();
    assert!(audit.transferred_mass_kg > 0.0);
    close(next.mass_kg(), expected.mass_kg(), 1e-17);
    close(next.energy_j(), expected.energy_j(), 1e-9);
    assert_eq!(next.exposed_area_m2(CellId(1)), Ok(0.0));
}

#[test]
fn newly_born_cell_has_no_premature_sublimation() {
    let thermal = scheduled(vec![birth(0.5, cell(1, 0.0, 0.0))], 100);
    let (next, _, audit) = thermal
        .advance_sublimating_to(
            CellId(1),
            &Ice,
            kinetics(),
            vacuum(),
            0.5,
            config(1.0),
            phase(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(next.mass_kg(), 0.02);
    assert_eq!(next.energy_j(), 0.0);
    assert_eq!(audit.phase_steps, 0);
    assert_eq!(audit.transferred_mass_kg, 0.0);
    assert_eq!(audit.thermal.incoming_mass_kg, 0.02);
}

#[test]
fn whole_advance_does_not_renew_the_thermal_work_budget() {
    let thermal = scheduled(vec![birth(0.0, cell(1, 0.0, 0.0))], 1)
        .advance_to(0.0, config(1.0), || Ok(()))
        .unwrap()
        .0;
    let original = thermal.cells().clone();
    let result = thermal.advance_sublimating_to(
        CellId(1),
        &Ice,
        kinetics(),
        vacuum(),
        1.0,
        config(0.25),
        phase(),
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(ThermalSublimationError::Thermal(DepositionError::Budget(_)))
    ));
    assert_eq!(thermal.cells(), &original);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn heat_phase_splitting_converges_under_time_refinement() {
    let thermal = initial(25.0);
    let run = |dt| {
        thermal
            .advance_sublimating_to(
                CellId(1),
                &Ice,
                kinetics(),
                vacuum(),
                4.0,
                config(dt),
                phase(),
                || Ok(()),
            )
            .unwrap()
            .0
            .cells()[&CellId(1)]
            .specific_h_j_kg
    };
    let reference = run(0.03125);
    let coarse = (run(1.0) - reference).abs();
    let medium = (run(0.5) - reference).abs();
    let fine = (run(0.25) - reference).abs();
    assert!(
        medium < coarse && fine < medium,
        "splitting errors: coarse={coarse:e}, medium={medium:e}, fine={fine:e}"
    );
}

#[test]
fn out_of_domain_trial_does_not_reject_a_converged_physical_endpoint() {
    let mut thermal = initial(0.0);
    let mut gas = vacuum();
    // The uncooled first trial predicts T<180 K, outside this caloric chart.
    // Latent feedback suppresses vapor pressure and admits a warmer endpoint.
    let step = thermal
        .sublimate_cell(
            CellId(1),
            &Ice,
            kinetics(),
            &mut gas,
            1000.0,
            phase(),
            || Ok(()),
        )
        .unwrap();
    assert!(step.surface_temperature_k > 180.0 && step.surface_temperature_k < 250.0);
    assert!(step.transferred_mass_kg > 0.0);
    close(thermal.energy_j() + step.vapor_energy_j, 0.0, 1e-7);
}

#[test]
fn empty_prebirth_time_jumps_with_zero_thermal_work_allowance() {
    let thermal = scheduled(vec![birth(1000.0, cell(1, 0.0, 0.0))], 0);
    let original_vapor = vacuum();
    let mut calls = 0;
    let (next, vapor, audit) = thermal
        .advance_sublimating_to(
            CellId(1),
            &Ice,
            kinetics(),
            original_vapor,
            1000.0,
            config(1e-12),
            phase(),
            || {
                calls += 1;
                if calls > 200 {
                    Err(DepositionError::Budget("unbounded prebirth interval loop"))
                } else {
                    Ok(())
                }
            },
        )
        .unwrap();
    assert_eq!(next.time_s(), 1000.0);
    assert_eq!(next.cells()[&CellId(1)].mass_kg, 0.02);
    assert_eq!(vapor, original_vapor);
    assert_eq!(audit.thermal.substeps, 0);
    assert_eq!(audit.thermal.events, 1);
    assert_eq!(audit.phase_steps, 0);
    assert_eq!(audit.transferred_mass_kg, 0.0);
}
