//! G0 shared inventory/atomicity and G1 actual multi-cell heat/phase coupling.
use super::*;
use crate::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, MaterialId, Reservoir, ThermalLaw,
    ThermalLink, TimedEvent,
};
use crate::sublimation::{SolverConfig, SublimationProperties};

struct Caloric(f64);
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 250.0 + h / self.0;
        if (180.0..=273.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.0
    }
}
static CALORIC: Caloric = Caloric(2000.0);
static SECOND_CALORIC: Caloric = Caloric(4000.0);

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
        let p = 100.0 * fs_math::det::exp(6000.0 * (1.0 / 250.0 - 1.0 / t));
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
fn selection() -> [ChamberSurface; 2] {
    [1, 2].map(|id| ChamberSurface {
        cell_id: CellId(id),
        kinetics: SurfaceKinetics {
            accommodation_coefficient: 0.01,
            mass_transfer_coefficient_m_s: None,
        },
    })
}
fn phase_config() -> ChamberSolverConfig {
    ChamberSolverConfig {
        surface_solver: SolverConfig {
            relative_tolerance: 1e-11,
            absolute_energy_tolerance_j: 1e-10,
            max_iterations: 96,
        },
        relative_mass_tolerance: 1e-10,
        absolute_mass_tolerance_kg: 1e-16,
        max_iterations: 96,
        max_surfaces: 4,
    }
}
fn thermal_config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-7,
    }
}
fn bath(pressure_pa: f64, volume_m3: f64) -> VaporReservoir {
    VaporReservoir::Finite {
        vapor_mass_kg: pressure_pa * volume_m3 * 0.018 / (8.314_462_618_153_24 * 250.0),
        volume_m3,
        temperature_k: 250.0,
    }
}
fn vapor_mass(bath: VaporReservoir) -> f64 {
    let VaporReservoir::Finite { vapor_mass_kg, .. } = bath else {
        panic!("finite chamber");
    };
    vapor_mass_kg
}
fn boundary(power: f64) -> Boundary {
    Boundary {
        ambient_k: 250.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: power,
    }
}
fn cell(id: u64, h: f64, power: f64) -> DepositedCell {
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
        boundary: boundary(power),
    }
}
fn birth(time_s: f64, cell: DepositedCell) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Birth(cell),
    }
}
fn scheduled(events: Vec<TimedEvent>, max_steps: usize) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &CALORIC), (MaterialId(2), &SECOND_CALORIC)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: max_steps,
        },
    )
    .unwrap()
}
fn initial(second_h: f64) -> DepositionSimulation<'static> {
    scheduled(
        vec![
            birth(0.0, cell(1, 0.0, 0.0)),
            birth(0.0, cell(2, second_h, 0.0)),
        ],
        10_000,
    )
    .advance_to(0.0, thermal_config(1.0), || Ok(()))
    .unwrap()
    .0
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn internal_energy_surface_refuses_the_entire_chamber_exchange() {
    struct InternalEnergy;
    impl ThermalLaw for InternalEnergy {
        fn storage_basis(&self) -> ThermalStorageBasis {
            ThermalStorageBasis::InternalEnergy
        }
        fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
            CALORIC.temperature_k(u)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            CALORIC.maximum_temperature_slope()
        }
    }
    let internal = InternalEnergy;
    let mut second = cell(2, -20_000.0, 0.0);
    second.material_id = MaterialId(2);
    let mut thermal = DepositionSimulation::new(
        vec![birth(0.0, cell(1, 0.0, 0.0)), birth(0.0, second)],
        &[(MaterialId(1), &CALORIC), (MaterialId(2), &internal)],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 4,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(1.0), || Ok(()))
    .unwrap()
    .0;
    let before = thermal.cells().clone();
    let before_mass = thermal.mass_kg();
    let before_energy = thermal.energy_j();
    let mut chamber = bath(40.0, 0.01);
    let before_chamber = chamber;
    let error = thermal
        .sublimate_chamber(
            &selection(),
            &Ice,
            &mut chamber,
            0.25,
            phase_config(),
            || Ok(()),
        )
        .unwrap_err();
    assert_eq!(
        error,
        ThermalSublimationError::Thermal(DepositionError::InvalidInput(
            "shared chamber sublimation requires enthalpy storage and an explicit open-mass model",
        )),
    );
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.mass_kg(), before_mass);
    assert_eq!(thermal.energy_j(), before_energy);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(chamber, before_chamber);
}

#[test]
fn identical_surfaces_fill_one_chamber_without_duplicate_inventory() {
    let mut thermal = initial(0.0);
    let old_mass = thermal.mass_kg();
    let old_energy = thermal.energy_j();
    let mut chamber = bath(0.0, 0.001);
    let step = thermal
        .sublimate_chamber(
            &selection(),
            &Ice,
            &mut chamber,
            5.0,
            phase_config(),
            || Ok(()),
        )
        .unwrap();
    close(thermal.mass_kg() + vapor_mass(chamber), old_mass, 1e-16);
    close(
        thermal.energy_j() + step.exchange.vapor_energy_j,
        old_energy,
        1e-7,
    );
    assert!(vapor_mass(chamber) > 0.0);
    assert!(vapor_mass(chamber) <= vapor_mass(bath(100.0, 0.001)));
    let a = step.exchange.surfaces[0];
    let b = step.exchange.surfaces[1];
    close(a.transferred_mass_kg, b.transferred_mass_kg, 1e-16);
    assert_eq!(a.reservoir, chamber);
    assert_eq!(b.reservoir, chamber);
    assert_eq!(step.cell_ids, vec![CellId(1), CellId(2)]);
}

#[test]
fn warm_solid_and_cold_solid_exchange_through_the_same_vapor() {
    let mut thermal = initial(-30_000.0);
    let old_energy = thermal.energy_j();
    let old_mass = thermal.mass_kg();
    let mut chamber = bath(50.0, 0.01);
    let old_vapor = vapor_mass(chamber);
    let step = thermal
        .sublimate_chamber(
            &selection(),
            &Ice,
            &mut chamber,
            0.25,
            phase_config(),
            || Ok(()),
        )
        .unwrap();
    assert!(step.exchange.surfaces[0].transferred_mass_kg > 0.0);
    assert!(step.exchange.surfaces[1].transferred_mass_kg < 0.0);
    assert!(thermal.temperature_k(CellId(1)).unwrap() < 250.0);
    assert!(thermal.temperature_k(CellId(2)).unwrap() > 235.0);
    close(
        thermal.mass_kg() + vapor_mass(chamber),
        old_mass + old_vapor,
        1e-16,
    );
    close(
        thermal.energy_j() - old_energy + step.exchange.vapor_energy_j,
        0.0,
        1e-7,
    );
}

#[test]
fn caller_order_cannot_change_cell_or_reservoir_results() {
    let mut left = initial(-20_000.0);
    let mut right = left.clone();
    let mut left_bath = bath(40.0, 0.01);
    let mut right_bath = left_bath;
    let forward = selection();
    let reverse = [forward[1], forward[0]];
    let a = left
        .sublimate_chamber(&forward, &Ice, &mut left_bath, 0.2, phase_config(), || {
            Ok(())
        })
        .unwrap();
    let b = right
        .sublimate_chamber(&reverse, &Ice, &mut right_bath, 0.2, phase_config(), || {
            Ok(())
        })
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(left.cells(), right.cells());
    assert_eq!(left_bath, right_bath);
}

#[test]
fn each_surface_uses_its_original_caloric_chart() {
    let mut second = cell(2, 0.0, 0.0);
    second.material_id = MaterialId(2);
    let mut thermal = scheduled(vec![birth(0.0, cell(1, 0.0, 0.0)), birth(0.0, second)], 100)
        .advance_to(0.0, thermal_config(1.0), || Ok(()))
        .unwrap()
        .0;
    thermal
        .sublimate_chamber(
            &selection(),
            &Ice,
            &mut bath(0.0, 0.01),
            0.25,
            phase_config(),
            || Ok(()),
        )
        .unwrap();
    assert!(thermal.temperature_k(CellId(1)).unwrap() < thermal.temperature_k(CellId(2)).unwrap());
    for (id, cp) in [(CellId(1), 2000.0), (CellId(2), 4000.0)] {
        close(
            thermal.temperature_k(id).unwrap(),
            250.0 + thermal.cells()[&id].specific_h_j_kg / cp,
            1e-12,
        );
    }
}

#[test]
fn original_heat_graph_and_joint_phase_close_the_complete_window() {
    let thermal = scheduled(
        vec![
            birth(0.0, cell(1, 0.0, 1.5)),
            birth(0.0, cell(2, -20_000.0, 1.5)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 7,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.0001,
                    conductance_w_k: 0.05,
                }),
            },
        ],
        10_000,
    )
    .advance_to(0.0, thermal_config(1.0), || Ok(()))
    .unwrap()
    .0;
    let old_bath = bath(40.0, 0.01);
    let (next, chamber, audit) = thermal
        .advance_chamber_sublimating_to(
            &selection(),
            &Ice,
            old_bath,
            0.2,
            thermal_config(0.025),
            phase_config(),
            || Ok(()),
        )
        .unwrap();
    close(audit.thermal.external_heat_j, 0.6, 1e-12);
    close(
        next.energy_j() - thermal.energy_j() + audit.vapor_energy_j,
        0.6,
        1e-7,
    );
    close(
        next.mass_kg() + vapor_mass(chamber),
        thermal.mass_kg() + vapor_mass(old_bath),
        1e-16,
    );
    close(
        audit
            .cells
            .values()
            .map(|row| row.transferred_mass_kg)
            .sum(),
        audit.transferred_mass_kg,
        1e-16,
    );
    close(
        audit.mass_balance_residual_kg,
        audit.phase_mass_residual_kg,
        1e-16,
    );
    assert!(audit.cells[&CellId(1)].recession_m > 0.0);
    assert!(audit.cells[&CellId(2)].recession_m < 0.0);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn final_checkpoint_refuses_before_any_of_the_cells_are_committed() {
    let original = initial(-20_000.0);
    let original_bath = bath(40.0, 0.01);
    let mut probe = original.clone();
    let mut probe_bath = original_bath;
    let mut checkpoints = 0;
    probe
        .sublimate_chamber(
            &selection(),
            &Ice,
            &mut probe_bath,
            0.2,
            phase_config(),
            || {
                checkpoints += 1;
                Ok(())
            },
        )
        .unwrap();
    let mut thermal = original.clone();
    let mut chamber = original_bath;
    let mut calls = 0;
    let result = thermal.sublimate_chamber(
        &selection(),
        &Ice,
        &mut chamber,
        0.2,
        phase_config(),
        || {
            calls += 1;
            if calls == checkpoints {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert_eq!(
        result,
        Err(ThermalSublimationError::Thermal(DepositionError::Cancelled))
    );
    assert_eq!(thermal.cells(), original.cells());
    assert_eq!(chamber, original_bath);
}

#[test]
fn endpoint_birth_is_not_exposed_before_it_exists() {
    let thermal = scheduled(
        vec![birth(0.0, cell(1, 0.0, 0.0)), birth(0.5, cell(2, 0.0, 0.0))],
        100,
    );
    let (next, _, audit) = thermal
        .advance_chamber_sublimating_to(
            &selection(),
            &Ice,
            bath(0.0, 0.001),
            0.5,
            thermal_config(1.0),
            phase_config(),
            || Ok(()),
        )
        .unwrap();
    assert!(audit.cells[&CellId(1)].transferred_mass_kg > 0.0);
    assert_eq!(audit.cells[&CellId(2)], CellPhaseAudit::default());
    assert_eq!(next.cells()[&CellId(2)].mass_kg, 0.02);
    assert_eq!(next.cells()[&CellId(2)].specific_h_j_kg, 0.0);
    assert_eq!(audit.last_step.unwrap().cell_ids, vec![CellId(1)]);
}

#[test]
fn endpoint_cover_does_not_rewrite_the_preceding_surface_exposure() {
    let mut covered = boundary(0.0);
    covered.reservoir = Some(Reservoir {
        area_m2: 0.001,
        conductance_w_k: 0.0,
        temperature_k: 250.0,
    });
    let thermal = scheduled(
        vec![
            birth(0.0, cell(1, 0.0, 0.0)),
            birth(0.0, cell(2, 0.0, 0.0)),
            TimedEvent {
                time_s: 0.5,
                event: Event::Boundary(CellId(1), covered),
            },
            TimedEvent {
                time_s: 0.5,
                event: Event::Boundary(CellId(2), covered),
            },
        ],
        100,
    )
    .advance_to(0.0, thermal_config(1.0), || Ok(()))
    .unwrap()
    .0;
    let mut expected = initial(0.0);
    let mut expected_bath = bath(0.0, 0.001);
    expected
        .sublimate_chamber(
            &selection(),
            &Ice,
            &mut expected_bath,
            0.5,
            phase_config(),
            || Ok(()),
        )
        .unwrap();
    let (next, chamber, _) = thermal
        .advance_chamber_sublimating_to(
            &selection(),
            &Ice,
            bath(0.0, 0.001),
            1.0,
            thermal_config(1.0),
            phase_config(),
            || Ok(()),
        )
        .unwrap();
    close(next.mass_kg(), expected.mass_kg(), 1e-16);
    close(next.energy_j(), expected.energy_j(), 1e-9);
    assert_eq!(chamber, expected_bath);
    assert_eq!(next.exposed_area_m2(CellId(1)), Ok(0.0));
    assert_eq!(next.exposed_area_m2(CellId(2)), Ok(0.0));
}

#[test]
fn duplicate_unknown_and_overbudget_selections_refuse_without_changes() {
    let mut thermal = initial(0.0);
    let before = thermal.cells().clone();
    let mut chamber = bath(0.0, 0.001);
    let original_bath = chamber;
    let a = selection()[0];
    for invalid in [
        vec![a, a],
        vec![ChamberSurface {
            cell_id: CellId(99),
            ..a
        }],
        Vec::new(),
    ] {
        assert!(
            thermal
                .sublimate_chamber(&invalid, &Ice, &mut chamber, 0.1, phase_config(), || Ok(()))
                .is_err()
        );
    }
    let mut restricted = phase_config();
    restricted.max_surfaces = 1;
    assert!(
        thermal
            .sublimate_chamber(&selection(), &Ice, &mut chamber, 0.1, restricted, || Ok(()))
            .is_err()
    );
    assert_eq!(thermal.cells(), &before);
    assert_eq!(chamber, original_bath);
}

#[test]
fn complete_advance_preserves_the_original_global_work_ceiling() {
    let thermal = scheduled(
        vec![birth(0.0, cell(1, 0.0, 0.0)), birth(0.0, cell(2, 0.0, 0.0))],
        1,
    )
    .advance_to(0.0, thermal_config(1.0), || Ok(()))
    .unwrap()
    .0;
    let before = thermal.cells().clone();
    let result = thermal.advance_chamber_sublimating_to(
        &selection(),
        &Ice,
        bath(0.0, 0.001),
        1.0,
        thermal_config(0.25),
        phase_config(),
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(ThermalSublimationError::Thermal(DepositionError::Budget(_)))
    ));
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn empty_chamber_timeline_jumps_to_births_without_phase_or_thermal_work() {
    let thermal = scheduled(
        vec![
            birth(1000.0, cell(1, 0.0, 0.0)),
            birth(1000.0, cell(2, 0.0, 0.0)),
        ],
        0,
    );
    let original_chamber = bath(0.0, 0.001);
    let mut calls = 0;
    let (next, chamber, audit) = thermal
        .advance_chamber_sublimating_to(
            &selection(),
            &Ice,
            original_chamber,
            1000.0,
            thermal_config(1e-12),
            phase_config(),
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
    assert_eq!(next.mass_kg(), 0.04);
    assert_eq!(chamber, original_chamber);
    assert_eq!(audit.thermal.substeps, 0);
    assert_eq!(audit.thermal.events, 2);
    assert_eq!(audit.phase_steps, 0);
    assert_eq!(audit.cells[&CellId(1)], CellPhaseAudit::default());
    assert_eq!(audit.cells[&CellId(2)], CellPhaseAudit::default());
}
