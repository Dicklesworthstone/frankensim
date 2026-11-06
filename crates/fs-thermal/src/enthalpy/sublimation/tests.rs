//! G0 conservation/transactionality and G1 coupling-refinement checks.
use super::*;
use crate::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, MaterialId, Reservoir,
    ThermalLaw, ThermalLink, ThermalStorageBasis, TimedEvent,
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
fn internal_energy_storage_refuses_sublimation_without_mutation() {
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
    let mut thermal = DepositionSimulation::new(
        vec![birth(0.0, cell(1, -20_000.0, 0.0))],
        &[(MaterialId(1), &internal)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 4,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0;
    let before = thermal.cells().clone();
    let before_mass = thermal.mass_kg();
    let before_energy = thermal.energy_j();
    let mut reservoir = VaporReservoir::Finite {
        vapor_mass_kg: 0.0,
        volume_m3: 0.001,
        temperature_k: 250.0,
    };
    let before_reservoir = reservoir;
    let error = thermal
        .sublimate_cell(
            CellId(1),
            &Ice,
            kinetics(),
            &mut reservoir,
            0.25,
            phase(),
            || Ok(()),
        )
        .unwrap_err();
    assert_eq!(
        error,
        ThermalSublimationError::Thermal(DepositionError::InvalidInput(
            "sublimation requires enthalpy storage and an explicit open-mass model",
        )),
    );
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.mass_kg(), before_mass);
    assert_eq!(thermal.energy_j(), before_energy);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(reservoir, before_reservoir);
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

#[test]
fn decimal_phase_intervals_do_not_create_a_roundoff_tail() {
    let thermal = scheduled(vec![birth(0.0, cell(1, 0.0, 0.0))], 32)
        .advance_to(0.0, config(1.0), || Ok(()))
        .unwrap()
        .0;
    let (next, _, audit) = thermal
        .advance_sublimating_to(
            CellId(1),
            &Ice,
            kinetics(),
            vacuum(),
            1.0,
            config(0.1),
            phase(),
            || Ok(()),
        )
        .unwrap();

    // Repeated addition stops first at 0.9999999999999999 and used to
    // invoke an eleventh phase solve for the remaining roundoff interval.
    assert_eq!(next.time_s(), 1.0);
    assert_eq!(audit.phase_steps, 10);
    assert!(audit.transferred_mass_kg > 0.0);
    assert!(next.temperature_k(CellId(1)).unwrap() < 250.0);
    close(
        next.mass_kg() - thermal.mass_kg() + audit.transferred_mass_kg,
        0.0,
        1e-16,
    );
    close(
        next.energy_j() - thermal.energy_j() + audit.vapor_energy_j,
        0.0,
        1e-7,
    );
    close(audit.energy_residual_j, 0.0, 1e-7);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn decimal_phase_epochs_restart_at_birth_and_cover_events() {
    let mut covered = boundary(0.0);
    covered.reservoir = Some(Reservoir {
        area_m2: 0.001,
        conductance_w_k: 0.0,
        temperature_k: 250.0,
    });
    let thermal = scheduled(
        vec![
            birth(0.35, cell(1, 0.0, 0.0)),
            TimedEvent {
                time_s: 0.85,
                event: Event::Boundary(CellId(1), covered),
            },
        ],
        32,
    );
    let run = |target| {
        thermal
            .advance_sublimating_to(
                CellId(1),
                &Ice,
                kinetics(),
                vacuum(),
                target,
                config(0.1),
                phase(),
                || Ok(()),
            )
            .unwrap()
    };
    let (at_cover, _, exposed) = run(0.85);
    let (next, _, audit) = run(1.35);

    // Empty time has no phase work. The five exposed intervals precede
    // coverage, and five covered intervals follow that new epoch.
    assert_eq!(exposed.phase_steps, 5);
    assert_eq!(audit.phase_steps, 10);
    assert_eq!(next.time_s(), 1.35);
    assert_eq!(audit.thermal.events, 2);
    assert!(exposed.transferred_mass_kg > 0.0);
    assert_eq!(next.cells(), at_cover.cells());
    assert_eq!(audit.transferred_mass_kg, exposed.transferred_mass_kg);
    assert_eq!(audit.latent_energy_j, exposed.latent_energy_j);
    assert_eq!(next.exposed_area_m2(CellId(1)), Ok(0.0));
    close(audit.energy_residual_j, 0.0, 1e-7);
}

#[test]
fn sublimation_and_deposition_balance_ignore_unchanged_neighbor_energy() {
    struct ShiftedCaloric(f64);
    impl ThermalLaw for ShiftedCaloric {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            CALORIC.temperature_k(h - self.0)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0 / 2000.0
        }
    }

    let neighbor_datum = 2.0_f64.powi(80);
    let neighbor_law = ShiftedCaloric(neighbor_datum);
    for birth_s in [0.0, 0.5] {
        for partial_pressure_pa in [0.0, 200.0] {
            let gas = VaporReservoir::Prescribed {
                partial_pressure_pa,
                temperature_k: 250.0,
            };
            let mut neighbor = cell(2, neighbor_datum, 0.0);
            neighbor.material_id = MaterialId(2);
            let thermal = DepositionSimulation::new(
                vec![
                    birth(0.0, neighbor.clone()),
                    birth(birth_s, cell(1, 0.0, 0.0)),
                ],
                &[
                    (MaterialId(1), &CALORIC as &dyn ThermalLaw),
                    (MaterialId(2), &neighbor_law as &dyn ThermalLaw),
                ],
                DepositionBudget {
                    max_cells: 2,
                    max_links: 0,
                    max_events: 2,
                    max_substeps: 16,
                },
            )
            .unwrap()
            .advance_to(0.0, config(1.0), || Ok(()))
            .unwrap()
            .0;
            let reference = scheduled(vec![birth(birth_s, cell(1, 0.0, 0.0))], 16)
                .advance_to(0.0, config(1.0), || Ok(()))
                .unwrap()
                .0;
            let (expected, expected_gas, expected_audit) = reference
                .advance_sublimating_to(
                    CellId(1),
                    &Ice,
                    kinetics(),
                    gas,
                    1.0,
                    config(0.125),
                    phase(),
                    || Ok(()),
                )
                .unwrap();
            let (next, next_gas, audit) = thermal
                .advance_sublimating_to(
                    CellId(1),
                    &Ice,
                    kinetics(),
                    gas,
                    1.0,
                    config(0.125),
                    phase(),
                    || Ok(()),
                )
                .unwrap();

            assert_eq!(next.cells()[&CellId(2)], neighbor);
            assert_eq!(next.cells()[&CellId(1)], expected.cells()[&CellId(1)]);
            assert_eq!(next_gas, expected_gas);
            assert_eq!(
                audit.transferred_mass_kg,
                expected_audit.transferred_mass_kg
            );
            assert_eq!(audit.latent_energy_j, expected_audit.latent_energy_j);
            assert_eq!(audit.vapor_energy_j, expected_audit.vapor_energy_j);
            assert_eq!(audit.phase_steps, expected_audit.phase_steps);
            assert_eq!(
                audit.thermal.incoming_mass_kg,
                expected_audit.thermal.incoming_mass_kg
            );
            close(audit.energy_residual_j, 0.0, 1e-7);
            close(audit.mass_balance_residual_kg, 0.0, 1e-16);
            close(
                audit.stored_energy_change_j,
                expected_audit.stored_energy_change_j,
                1e-12,
            );
            assert_eq!(audit.transferred_mass_kg > 0.0, partial_pressure_pa == 0.0);
            assert_eq!(next.energy_j(), thermal.energy_j());
            assert!(audit.vapor_energy_j.abs() > 0.01);
            assert!(audit.stored_energy_change_j.abs() > 0.01);
            assert_eq!(thermal.time_s(), 0.0);
        }
    }
}

#[test]
fn supplied_heat_drives_native_phase_transfer_without_sensible_overshoot() {
    let duration = 1000.0;
    for partial_pressure_pa in [0.0, 200.0] {
        let mut thermal = initial(7.0);
        let original = thermal.cells()[&CellId(1)].clone();
        let mut gas = VaporReservoir::Prescribed {
            partial_pressure_pa,
            temperature_k: 250.0,
        };
        // Independent Hertz-Knudsen flux at the fixture's original 250 K.
        let transfer = (0.018
            / (2.0 * core::f64::consts::PI * crate::sublimation::MOLAR_GAS_CONSTANT * 250.0))
            .sqrt()
            * 0.01
            * (100.0 - partial_pressure_pa)
            * original.geometry.surface_area_m2
            * duration;
        let latent = crate::sublimation::MOLAR_GAS_CONSTANT / 0.018 * 6000.0;
        let heat = transfer * latent;
        let sensible_only_t = 250.0 + heat / (original.mass_kg * 2000.0);
        assert!(!(180.0..=273.0).contains(&sensible_only_t));
        let step = thermal
            .sublimate_cell_with_heat(
                CellId(1),
                &Ice,
                kinetics(),
                &mut gas,
                duration,
                heat,
                phase(),
                || Ok(()),
            )
            .unwrap();

        assert_eq!(thermal.temperature_k(CellId(1)).unwrap(), 250.0);
        assert_eq!(step.supplied_heat_j, heat);
        close(step.transferred_mass_kg, transfer, 1e-16);
        close(step.latent_energy_j, heat, 1e-7);
        close(
            thermal.mass_kg() + step.transferred_mass_kg,
            original.mass_kg,
            1e-17,
        );
        close(thermal.energy_j() + step.vapor_energy_j, heat, 1e-7);
        assert_eq!(step.transferred_mass_kg > 0.0, partial_pressure_pa == 0.0);
        assert_eq!(thermal.cells()[&CellId(1)].geometry, original.geometry);
        assert_eq!(thermal.cells()[&CellId(1)].boundary, original.boundary);
        assert_eq!(thermal.time_s(), 0.0);
        // The authored 7 W is retained for the heat owner. This direct operation
        // spends exactly its supplied heat, without applying that source again.
    }
}

#[test]
fn heated_sublimation_refuses_invalid_heat_and_final_cancellation_atomically() {
    let reservoir = VaporReservoir::Finite {
        vapor_mass_kg: 0.0,
        volume_m3: 0.1,
        temperature_k: 250.0,
    };
    let mut successful = initial(0.0);
    let mut successful_gas = reservoir;
    let mut polls = 0;
    successful
        .sublimate_cell_with_heat(
            CellId(1),
            &Ice,
            kinetics(),
            &mut successful_gas,
            1.0,
            10.0,
            phase(),
            || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();

    for cancel_at in [1, polls] {
        let mut thermal = initial(0.0);
        let before = thermal.cells().clone();
        let mut gas = reservoir;
        let mut visited = 0;
        let result = thermal.sublimate_cell_with_heat(
            CellId(1),
            &Ice,
            kinetics(),
            &mut gas,
            1.0,
            10.0,
            phase(),
            || {
                visited += 1;
                if visited == cancel_at {
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
        assert_eq!(thermal.cells(), &before);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(gas, reservoir);
    }

    let mut thermal = initial(0.0);
    let before = thermal.cells().clone();
    for heat in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut gas = reservoir;
        assert!(matches!(
            thermal.sublimate_cell_with_heat(
                CellId(1),
                &Ice,
                kinetics(),
                &mut gas,
                1.0,
                heat,
                phase(),
                || Ok(()),
            ),
            Err(ThermalSublimationError::Surface(SublimationError::Invalid(
                _
            )))
        ));
        assert_eq!(thermal.cells(), &before);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(gas, reservoir);
    }
}

#[test]
fn coupled_authored_heat_sustains_latent_transfer_without_sensible_overshoot() {
    let duration = 1000.0;
    for pressure in [0.0, 200.0] {
        let transfer = (0.018
            / (2.0 * core::f64::consts::PI * crate::sublimation::MOLAR_GAS_CONSTANT * 250.0))
            .sqrt()
            * 0.01
            * (100.0 - pressure)
            * 0.001
            * duration;
        let heat = transfer * crate::sublimation::MOLAR_GAS_CONSTANT / 0.018 * 6000.0;
        let thermal = initial(heat / duration);
        let before = thermal.cells().clone();
        let gas = VaporReservoir::Prescribed {
            partial_pressure_pa: pressure,
            temperature_k: 250.0,
        };
        // A full heat-only intermediate lies outside the admitted solid chart.
        assert!(
            thermal
                .advance_sublimating_to(
                    CellId(1),
                    &Ice,
                    kinetics(),
                    gas,
                    duration,
                    config(duration),
                    phase(),
                    || Ok(()),
                )
                .is_err()
        );
        let (next, returned, audit) = thermal
            .advance_coupled_sublimating_to(
                CellId(1),
                &Ice,
                kinetics(),
                gas,
                duration,
                config(duration),
                phase(),
                || Ok(()),
            )
            .unwrap();
        assert_eq!(next.temperature_k(CellId(1)), Ok(250.0));
        assert_eq!(next.time_s(), duration);
        assert_eq!(audit.substeps, 1);
        assert_eq!(audit.phase_steps, 1);
        assert_eq!(audit.events, 0);
        assert_eq!(returned, gas);
        close(audit.transferred_mass_kg, transfer, 1e-16);
        close(audit.external_heat_j, heat, 1e-10);
        close(audit.latent_energy_j, heat, 1e-8);
        close(audit.last_step.unwrap().supplied_heat_j, heat, 1e-10);
        close(audit.energy_residual_j, 0.0, 1e-8);
        close(audit.mass_balance_residual_kg, 0.0, 1e-17);
        close(
            next.energy_j() - thermal.energy_j() + audit.vapor_energy_j,
            heat,
            1e-8,
        );
        assert_eq!(
            next.cells()[&CellId(1)].geometry,
            before[&CellId(1)].geometry
        );
        assert_eq!(
            next.cells()[&CellId(1)].boundary,
            before[&CellId(1)].boundary
        );
        assert_eq!(thermal.cells(), &before);
        assert_eq!(thermal.time_s(), 0.0);
    }
}

#[test]
fn coupled_source_baths_and_contact_debit_one_shared_heat_interval() {
    let mut surface = cell(1, 0.0, 2.0);
    surface.boundary.ambient_k = 260.0;
    surface.boundary.htc_w_m2_k = 10.0;
    surface.boundary.reservoir = Some(Reservoir {
        area_m2: 0.0002,
        conductance_w_k: 0.2,
        temperature_k: 270.0,
    });
    let thermal = scheduled(
        vec![
            birth(0.0, surface),
            birth(0.0, cell(2, 20_000.0, 0.0)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 9,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.0001,
                    conductance_w_k: 0.3,
                }),
            },
        ],
        4,
    )
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0;
    let (next, _, audit) = thermal
        .advance_coupled_sublimating_to(
            CellId(1),
            &Ice,
            kinetics(),
            vacuum(),
            0.25,
            config(0.25),
            phase(),
            || Ok(()),
        )
        .unwrap();
    let step = audit.last_step.unwrap();
    // Start-state powers: source 2, convection .07, bath 4, contact 3 W.
    close(step.supplied_heat_j, 2.2675, 1e-14);
    close(audit.external_heat_j, 1.5175, 1e-14);
    close(
        next.cells()[&CellId(2)].specific_h_j_kg,
        20_000.0 - 0.75 / 0.02,
        1e-12,
    );
    close(step.solid.area_m2, 0.0007, 1e-18);
    assert_eq!(audit.substeps, 1);
    close(audit.energy_residual_j, 0.0, 1e-7);
    close(
        next.energy_j() - thermal.energy_j() + audit.vapor_energy_j,
        1.5175,
        1e-7,
    );

    // Independent pressure/kinetic expression and scalar energy balance.
    let residual = |t: f64| {
        let p = 100.0 * (6000.0 * (1.0 / 250.0 - 1.0 / t)).exp();
        let dm = 0.01
            * (0.018 / (2.0 * core::f64::consts::PI * 8.314_462_618_153_24 * t)).sqrt()
            * p
            * 0.0007
            * 0.25;
        (
            0.02 * 2000.0 * (t - 250.0) + dm * 8.314_462_618_153_24 / 0.018 * 6000.0 - 2.2675,
            dm,
        )
    };
    let (mut lower, mut upper) = (250.0, 251.0);
    for _ in 0..64 {
        let midpoint = lower + (upper - lower) * 0.5;
        if residual(midpoint).0 > 0.0 {
            upper = midpoint;
        } else {
            lower = midpoint;
        }
    }
    let expected_t = lower + (upper - lower) * 0.5;
    close(step.surface_temperature_k, expected_t, 1e-9);
    close(step.transferred_mass_kg, residual(expected_t).1, 1e-14);
}

#[test]
fn coupled_birth_coverage_and_global_budget_remain_atomic() {
    let flux = 0.01
        * (0.018 / (2.0 * core::f64::consts::PI * 8.314_462_618_153_24 * 250.0_f64)).sqrt()
        * 100.0
        * 0.001;
    let power = flux * 8.314_462_618_153_24 / 0.018 * 6000.0;
    let mut covered = boundary(0.0);
    covered.reservoir = Some(Reservoir {
        area_m2: 0.001,
        conductance_w_k: 0.0,
        temperature_k: 250.0,
    });
    let events = vec![
        birth(0.25, cell(1, 0.0, power)),
        TimedEvent {
            time_s: 0.75,
            event: Event::Boundary(CellId(1), covered),
        },
    ];
    let thermal = scheduled(events.clone(), 6);
    let mut polls = 0;
    let (next, _, audit) = thermal
        .advance_coupled_sublimating_to(
            CellId(1),
            &Ice,
            kinetics(),
            vacuum(),
            1.0,
            config(0.125),
            phase(),
            || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();
    assert_eq!(next.time_s(), 1.0);
    assert_eq!(audit.substeps, 6);
    assert_eq!(audit.phase_steps, 6);
    assert_eq!(audit.events, 2);
    assert_eq!(next.exposed_area_m2(CellId(1)), Ok(0.0));
    assert_eq!(next.cells()[&CellId(1)].boundary.source_w, 0.0);
    close(audit.external_heat_j, power * 0.5, 1e-13);
    close(audit.transferred_mass_kg, flux * 0.5, 1e-17);
    close(audit.energy_residual_j, 0.0, 1e-8);

    let (at_birth, _, birth_audit) = scheduled(events.clone(), 0)
        .advance_coupled_sublimating_to(
            CellId(1),
            &Ice,
            kinetics(),
            vacuum(),
            0.25,
            config(0.125),
            phase(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(at_birth.cells().len(), 1);
    assert_eq!(birth_audit.substeps, 0);
    assert_eq!(birth_audit.phase_steps, 0);
    assert_eq!(birth_audit.incoming_mass_kg, 0.02);

    let limited = scheduled(events, 5);
    assert!(matches!(
        limited.advance_coupled_sublimating_to(
            CellId(1),
            &Ice,
            kinetics(),
            vacuum(),
            1.0,
            config(0.125),
            phase(),
            || Ok(()),
        ),
        Err(ThermalSublimationError::Thermal(DepositionError::Budget(_)))
    ));
    assert!(limited.cells().is_empty());
    assert_eq!(limited.time_s(), 0.0);

    let mut visits = 0;
    assert!(matches!(
        thermal.advance_coupled_sublimating_to(
            CellId(1),
            &Ice,
            kinetics(),
            vacuum(),
            1.0,
            config(0.125),
            phase(),
            || {
                visits += 1;
                if visits == polls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        ),
        Err(ThermalSublimationError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(visits, polls);
    assert!(thermal.cells().is_empty());
    assert_eq!(thermal.time_s(), 0.0);
}
