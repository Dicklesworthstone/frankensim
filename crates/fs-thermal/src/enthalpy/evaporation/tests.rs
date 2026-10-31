//! G0 conservation/rollback and G1 complete heat/evaporation refinement.
use super::*;
use crate::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, MaterialId, Reservoir, ThermalLaw,
    ThermalLink, TimedEvent,
};
use crate::evaporation::EvaporationProperties;

const CP: f64 = 4180.0;
const LATENT: f64 = 2.4e6;
struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h / CP;
        if (280.0..=330.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / CP
    }
}
static CALORIC: Caloric = Caloric;

// Analytic liquid fixture; not an empirical water accuracy claim.
struct Liquid;
impl EvaporationLaw for Liquid {
    fn temperature_range_k(&self) -> [f64; 2] {
        [280.0, 330.0]
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        if !(280.0..=330.0).contains(&t) {
            return Err(EvaporationError::TemperatureOutOfRange {
                temperature_k: t,
                minimum_k: 280.0,
                maximum_k: 330.0,
            });
        }
        let beta = LATENT * 0.018 / 8.314_462_618_153_24;
        let p = 3600.0 * fs_math::det::exp(beta * (1.0 / 300.0 - 1.0 / t));
        Ok(EvaporationProperties {
            saturation_pressure_pa: p,
            saturation_pressure_derivative_pa_k: p * beta / (t * t),
            latent_heat_j_kg: LATENT,
            liquid_density_kg_m3: 1000.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}
fn film(p: f64) -> GasFilm {
    GasFilm {
        total_pressure_pa: 100_000.0,
        vapor_partial_pressure_pa: p,
        temperature_k: 300.0,
        mass_transfer_coefficient_m_s: 0.01,
    }
}
fn phase() -> SolverConfig {
    SolverConfig {
        specific_enthalpy_bracket_j_kg: [-20.0 * CP, 30.0 * CP],
        relative_tolerance: 1e-12,
        absolute_energy_tolerance_j: 1e-9,
        max_iterations: 128,
    }
}
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-6,
    }
}
fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
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
            volume_m3: 2e-5,
            surface_area_m2: 0.01,
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
fn initial(source_w: f64) -> DepositionSimulation<'static> {
    scheduled(vec![birth(0.0, cell(1, 0.0, source_w))], 10_000)
        .advance_to(0.0, config(1.0), || Ok(()))
        .unwrap()
        .0
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn evaporation_cools_and_condensation_warms_actual_liquid_inventory() {
    for (pressure, direction) in [(0.0, 1.0), (6000.0, -1.0)] {
        let mut thermal = initial(0.0);
        let original = thermal.cells()[&CellId(1)].clone();
        let step = thermal
            .evaporate_cell(CellId(1), &Liquid, film(pressure), 1.0, phase(), || Ok(()))
            .unwrap();
        assert!(direction * step.transferred_mass_kg > 0.0);
        assert!(direction * (300.0 - thermal.temperature_k(CellId(1)).unwrap()) > 0.0);
        close(
            thermal.mass_kg() + step.transferred_mass_kg,
            original.mass_kg,
            1e-17,
        );
        close(thermal.energy_j() + step.vapor_energy_j, 0.0, 1e-8);
        close(
            step.latent_energy_j,
            step.transferred_mass_kg * LATENT,
            1e-10,
        );
        assert_eq!(thermal.cells()[&CellId(1)].geometry, original.geometry);
        assert_eq!(thermal.cells()[&CellId(1)].boundary, original.boundary);
        assert_eq!(thermal.time_s(), 0.0);
    }
    // Ten decimal-sized physical intervals must not leave a tiny extra
    // evaporation step whose nonzero mass cannot be represented.
    let (endpoint, audit) = initial(0.0)
        .advance_evaporating_to(
            CellId(1),
            &Liquid,
            film(0.0),
            2.0,
            config(0.2),
            phase(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(endpoint.time_s(), 2.0);
    assert_eq!(audit.phase_steps, 10);
}

#[test]
fn conduction_source_and_evaporation_close_one_whole_window_account() {
    let thermal = scheduled(
        vec![
            birth(0.0, cell(1, 0.0, 2.0)),
            birth(0.0, cell(2, 10.0 * CP, 0.0)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 1,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.001,
                    conductance_w_k: 0.3,
                }),
            },
        ],
        1000,
    )
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0;
    let (next, audit) = thermal
        .advance_evaporating_to(
            CellId(1),
            &Liquid,
            film(0.0),
            4.0,
            config(0.25),
            phase(),
            || Ok(()),
        )
        .unwrap();
    assert!(audit.transferred_mass_kg > 0.0);
    assert!(next.temperature_k(CellId(2)).unwrap() < 310.0);
    close(audit.thermal.external_heat_j, 8.0, 1e-12);
    close(
        next.mass_kg() - thermal.mass_kg() + audit.transferred_mass_kg,
        0.0,
        1e-16,
    );
    close(
        next.energy_j() - thermal.energy_j() - 8.0 + audit.vapor_energy_j,
        0.0,
        1e-6,
    );
    assert!(audit.energy_residual_j.abs() < 1e-6);
    assert_eq!(audit.phase_steps, 16);
    assert_eq!(next.links(), thermal.links());
}

#[test]
fn endpoint_covering_keeps_start_exposure_and_birth_has_no_premature_evaporation() {
    let mut covered = boundary(0.0);
    covered.reservoir = Some(Reservoir {
        area_m2: 0.01,
        conductance_w_k: 0.0,
        temperature_k: 300.0,
    });
    let thermal = scheduled(
        vec![
            birth(0.5, cell(1, 0.0, 0.0)),
            TimedEvent {
                time_s: 1.0,
                event: Event::Boundary(CellId(1), covered),
            },
        ],
        100,
    );
    let (born, audit) = thermal
        .advance_evaporating_to(
            CellId(1),
            &Liquid,
            film(0.0),
            0.5,
            config(1.0),
            phase(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(born.mass_kg(), 0.02);
    assert_eq!(audit.phase_steps, 0);
    let mut expected = initial(0.0);
    expected
        .evaporate_cell(CellId(1), &Liquid, film(0.0), 0.5, phase(), || Ok(()))
        .unwrap();
    let (next, audit) = thermal
        .advance_evaporating_to(
            CellId(1),
            &Liquid,
            film(0.0),
            2.0,
            config(1.0),
            phase(),
            || Ok(()),
        )
        .unwrap();
    close(next.mass_kg(), expected.mass_kg(), 1e-17);
    close(next.energy_j(), expected.energy_j(), 1e-8);
    assert_eq!(next.exposed_area_m2(CellId(1)), Ok(0.0));
    assert!(audit.transferred_mass_kg > 0.0);
    assert_eq!(audit.last_step.unwrap().transferred_mass_kg, 0.0);
}

#[test]
fn pure_liquid_branch_and_total_work_budget_refuse_atomically() {
    let thermal = initial(1000.0);
    let before = thermal.cells().clone();
    let mut narrow = phase();
    narrow.specific_enthalpy_bracket_j_kg = [-10_000.0, 1000.0];
    assert!(matches!(
        thermal.advance_evaporating_to(
            CellId(1),
            &Liquid,
            film(0.0),
            1.0,
            config(1.0),
            narrow,
            || Ok(())
        ),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(thermal.cells(), &before);
    let thermal = scheduled(vec![birth(0.0, cell(1, 0.0, 0.0))], 1)
        .advance_to(0.0, config(1.0), || Ok(()))
        .unwrap()
        .0;
    let before = thermal.cells().clone();
    assert!(matches!(
        thermal.advance_evaporating_to(
            CellId(1),
            &Liquid,
            film(0.0),
            1.0,
            config(0.25),
            phase(),
            || Ok(())
        ),
        Err(ThermalEvaporationError::Thermal(DepositionError::Budget(_)))
    ));
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn late_cancellation_discards_the_complete_heat_and_phase_candidate() {
    let thermal = initial(2.0);
    let mut checkpoints = 0;
    thermal
        .advance_evaporating_to(
            CellId(1),
            &Liquid,
            film(0.0),
            1.0,
            config(0.25),
            phase(),
            || {
                checkpoints += 1;
                Ok(())
            },
        )
        .unwrap();
    let mut seen = 0;
    let before = thermal.cells().clone();
    let result = thermal.advance_evaporating_to(
        CellId(1),
        &Liquid,
        film(0.0),
        1.0,
        config(0.25),
        phase(),
        || {
            seen += 1;
            if seen == checkpoints {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        result,
        Err(ThermalEvaporationError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn invalid_film_and_closed_internal_energy_storage_refuse_without_mutation() {
    let mut thermal = initial(0.0);
    let before = thermal.cells().clone();
    assert!(
        thermal
            .evaporate_cell(CellId(1), &Liquid, film(100_001.0), 1.0, phase(), || Ok(()))
            .is_err()
    );
    assert_eq!(thermal.cells(), &before);
    // A prescribed dry bath cannot erase the cell when cooling leaves the
    // admitted liquid domain or exhausts its finite inventory.
    assert!(
        thermal
            .evaporate_cell(CellId(1), &Liquid, film(0.0), 1e9, phase(), || Ok(()))
            .is_err()
    );
    assert_eq!(thermal.cells(), &before);
    struct Internal;
    impl ThermalLaw for Internal {
        fn storage_basis(&self) -> ThermalStorageBasis {
            ThermalStorageBasis::InternalEnergy
        }
        fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
            CALORIC.temperature_k(u)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0 / CP
        }
    }
    let law = Internal;
    let mut thermal = DepositionSimulation::new(
        vec![birth(0.0, cell(1, 0.0, 0.0))],
        &[(MaterialId(1), &law)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 10,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0;
    let before = thermal.cells().clone();
    assert!(matches!(
        thermal.evaporate_cell(CellId(1), &Liquid, film(0.0), 1.0, phase(), || Ok(())),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput(_)
        ))
    ));
    assert_eq!(thermal.cells(), &before);
}

#[test]
fn complete_heat_phase_split_converges_under_time_refinement() {
    let thermal = initial(2.0);
    let run = |dt| {
        thermal
            .advance_evaporating_to(
                CellId(1),
                &Liquid,
                film(0.0),
                8.0,
                config(dt),
                phase(),
                || Ok(()),
            )
            .unwrap()
            .0
            .temperature_k(CellId(1))
            .unwrap()
    };
    let reference = run(0.015625);
    let coarse = (run(1.0) - reference).abs();
    let medium = (run(0.5) - reference).abs();
    let fine = (run(0.25) - reference).abs();
    assert!(
        coarse > 1e-6,
        "refinement must resolve nonzero coupling error"
    );
    assert!(medium < 0.7 * coarse, "{coarse:e}, {medium:e}");
    assert!(fine < 0.7 * medium, "{medium:e}, {fine:e}");
}

fn external(thermal: &DepositionSimulation<'_>, id: u64) -> LiquidCellExchange {
    LiquidCellExchange {
        expected: thermal.cells()[&CellId(id)].clone(),
        interval_exposed_area_m2: 0.01,
        replacement: LiquidSurfaceState {
            mass_kg: 0.0199,
            specific_enthalpy_j_kg: -1000.0,
            area_m2: 0.01,
        },
        pure_liquid_enthalpy_range_j_kg: phase().specific_enthalpy_bracket_j_kg,
    }
}

#[test]
fn externally_solved_liquid_endpoint_preserves_original_owner_and_rejects_staleness() {
    let mut thermal = initial(2.0);
    let update = external(&thermal, 1);
    let mut expected = update.expected.clone();
    expected.mass_kg = update.replacement.mass_kg;
    expected.specific_h_j_kg = update.replacement.specific_enthalpy_j_kg;
    thermal
        .apply_liquid_exchanges(&[update.clone()], || Ok(()))
        .unwrap();
    assert_eq!(thermal.cells()[&CellId(1)], expected);
    assert_eq!(thermal.time_s(), 0.0);
    let accepted = thermal.cells().clone();
    assert!(
        thermal
            .apply_liquid_exchanges(&[update], || Ok(()))
            .is_err()
    );
    assert_eq!(thermal.cells(), &accepted);
}

#[test]
fn external_liquid_batch_validates_all_cells_and_final_cancellation_before_writing() {
    let mut thermal = scheduled(
        vec![birth(0.0, cell(1, 0.0, 0.0)), birth(0.0, cell(2, 0.0, 0.0))],
        10,
    )
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0;
    let before = thermal.cells().clone();
    let updates = [external(&thermal, 1), external(&thermal, 2)];
    for which in 0..4 {
        let mut bad = updates.clone();
        match which {
            0 => bad[1].replacement.mass_kg = 0.0,
            1 => bad[1].replacement.specific_enthalpy_j_kg = 1e7,
            2 => bad[1].expected.geometry.surface_area_m2 *= 2.0,
            _ => bad[1] = bad[0].clone(),
        }
        assert!(thermal.apply_liquid_exchanges(&bad, || Ok(())).is_err());
        assert_eq!(thermal.cells(), &before);
    }
    let mut polls = 0;
    assert_eq!(
        thermal.apply_liquid_exchanges(&updates, || {
            polls += 1;
            if polls == 3 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    );
    assert_eq!(thermal.cells(), &before);
}
