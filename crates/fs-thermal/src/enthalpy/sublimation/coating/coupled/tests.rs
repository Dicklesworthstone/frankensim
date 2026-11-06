//! G1/G4: authored latent heat, one native heat graph, and atomic event history.
use super::super::{ChamberCoatingSurface, CoupledCoatingAudit};
use super::*;
use crate::sublimation::chamber::ChamberSolverConfig;

fn chamber_controls() -> ChamberSolverConfig {
    ChamberSolverConfig {
        surface_solver: phase(),
        relative_mass_tolerance: 1e-11,
        absolute_mass_tolerance_kg: 1e-16,
        max_iterations: 128,
        max_surfaces: 4,
    }
}

fn finite_bath(pressure_pa: f64) -> VaporReservoir {
    VaporReservoir::Finite {
        vapor_mass_kg: pressure_pa * 0.001 * 0.018 / (MOLAR_GAS_CONSTANT * INITIAL_T),
        volume_m3: 0.001,
        temperature_k: INITIAL_T,
    }
}

fn gas_mass(reservoir: VaporReservoir) -> f64 {
    match reservoir {
        VaporReservoir::Finite { vapor_mass_kg, .. } => vapor_mass_kg,
        _ => panic!("expected the single finite chamber"),
    }
}

#[test]
fn authored_heating_and_cooling_feed_latent_transfer_without_an_intermediate_domain_exit() {
    let duration = 1000.0;
    for pressure in [0.0, 20.0] {
        let rate = if pressure == 0.0 {
            growth_per_second()
        } else {
            -growth_per_second()
        };
        let initial_deposit = if pressure == 0.0 { 0.02 } else { 0.0 };
        let source = rate * LATENT;
        let initial = simulation(
            vec![birth(0.0, cell(1, INITIAL_T, source))],
            &CALORIC,
            &[surface(1, initial_deposit, 1.0, bath(pressure), &deposit_h)],
            1,
        );
        let before = initial.clone();
        let sensible_only = INITIAL_T
            + source * duration / (SUBSTRATE_MASS * SUBSTRATE_CP + initial_deposit * DEPOSIT_CP);
        assert!(!(100.0..=400.0).contains(&sensible_only));
        assert!(
            initial
                .advance_to(duration, config(duration), phase(), || Ok(()))
                .is_err()
        );

        let (next, audit) = initial
            .advance_coupled_to(duration, config(duration), phase(), || Ok(()))
            .unwrap();
        let transfer = rate * duration;
        close(next.temperature_k(CellId(1)).unwrap(), INITIAL_T, 1e-9);
        close(
            next.coatings()[&CellId(1)].mass_kg,
            initial_deposit - transfer,
            1e-16,
        );
        close(audit.transferred_mass_kg, transfer, 1e-16);
        close(audit.latent_energy_j, source * duration, 1e-7);
        close(audit.external_heat_j, source * duration, 1e-7);
        close(
            audit.last_steps[&CellId(1)].phase.supplied_heat_j,
            source * duration,
            1e-7,
        );
        close(
            audit.stored_energy_change_j + audit.vapor_energy_j - audit.external_heat_j,
            0.0,
            1e-6,
        );
        close(audit.energy_residual_j, 0.0, 1e-6);
        assert!(audit.absolute_energy_residual_j <= 1e-6);
        assert_eq!(audit.substeps, 1);
        assert_eq!(audit.phase_steps, 1);
        assert_eq!(next.cells()[&CellId(1)].mass_kg, SUBSTRATE_MASS);
        assert_eq!(
            next.cells()[&CellId(1)].geometry,
            initial.cells()[&CellId(1)].geometry
        );
        assert_eq!(next.time_s(), duration);
        same_state(&initial, &before);
    }
}

#[test]
fn shared_heat_keeps_native_contact_debit_and_an_unheated_bare_receiver() {
    // Unlike the constant-flux fixture, this source permits warm material to
    // supply a colder bare surface through their one common vapor pressure.
    struct TemperatureDependentVapor;
    impl SublimationLaw for TemperatureDependentVapor {
        fn temperature_range_k(&self) -> [f64; 2] {
            [100.0, 400.0]
        }
        fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
            let p = 10.0 * (0.06 * (t - INITIAL_T)).exp();
            Ok(SublimationProperties {
                saturation_pressure_pa: p,
                saturation_pressure_derivative_pa_k: 0.06 * p,
                latent_heat_j_kg: LATENT,
                latent_heat_derivative_j_kg_k: 0.0,
                solid_density_kg_m3: 900.0,
                molar_mass_kg_mol: 0.018,
            })
        }
    }
    let law = TemperatureDependentVapor;
    let thermal = DepositionSimulation::new(
        vec![
            birth(0.0, cell(1, 250.0, 10.0)),
            birth(0.0, cell(2, 220.0, 0.0)),
            birth(0.0, cell(3, 240.0, 0.0)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 8,
                    a: CellId(1),
                    b: CellId(3),
                    area_m2: 0.0002,
                    conductance_w_k: 0.1,
                }),
            },
        ],
        &[(MaterialId(1), &CALORIC)],
        DepositionBudget {
            max_cells: 3,
            max_links: 1,
            max_events: 4,
            max_substeps: 1,
        },
    )
    .unwrap();
    let kinetics = SurfaceKinetics {
        accommodation_coefficient: 1.0,
        mass_transfer_coefficient_m_s: None,
    };
    let reservoir = finite_bath(10.0);
    let initial = CoatingSimulation::new_shared_chamber(
        thermal,
        &[
            ChamberCoatingSurface {
                cell_id: CellId(2),
                kinetics,
                initial_mass_kg: 0.0,
            },
            ChamberCoatingSurface {
                cell_id: CellId(1),
                kinetics,
                initial_mass_kg: 0.001,
            },
        ],
        &law,
        &deposit_h,
        reservoir,
        || Ok(()),
    )
    .unwrap();
    let (next, audit) = initial
        .advance_coupled_chamber_to(0.25, config(0.25), chamber_controls(), || Ok(()))
        .unwrap();
    let donor = audit.last_steps[&CellId(1)];
    let receiver = audit.last_steps[&CellId(2)];
    // Native contact sends 1 W from the hot donor to the ordinary finite body.
    // The source supplies 10 W: 2.25 J reaches phase and 0.25 J warms that body.
    close(donor.phase.supplied_heat_j, 2.25, 1e-14);
    assert_eq!(receiver.phase.supplied_heat_j, 0.0);
    close(audit.external_heat_j, 2.5, 1e-14);
    close(
        next.cells()[&CellId(3)].specific_h_j_kg,
        240.0 * SUBSTRATE_CP + 0.25 / SUBSTRATE_MASS,
        1e-10,
    );
    close(donor.state.area_m2, 0.0008, 1e-18);
    assert!(donor.phase.transferred_mass_kg > 0.0);
    assert!(receiver.phase.transferred_mass_kg < 0.0);
    assert!(next.coatings()[&CellId(2)].mass_kg > 0.0);
    let final_gas = next.shared_reservoir().unwrap();
    for coating in next.coatings().values() {
        assert_eq!(coating.reservoir, final_gas);
    }
    close(
        next.deposit_mass_kg() + gas_mass(final_gas),
        0.001 + gas_mass(reservoir),
        1e-14,
    );
    close(
        audit.stored_energy_change_j - audit.incoming_energy_j + audit.vapor_energy_j
            - audit.external_heat_j,
        0.0,
        1e-6,
    );
    close(audit.energy_residual_j, 0.0, 1e-6);
    assert!(audit.absolute_energy_residual_j <= 1e-6);
    assert_eq!(audit.substeps, 1);
    assert_eq!(audit.phase_steps, 2);
    assert_eq!(audit.events, 4);
    assert!(audit.pressure_iterations > 0);
    assert_eq!(next.time_s(), 0.25);
}

fn event_history(shared: bool) -> CoatingSimulation<'static> {
    let mut covered = boundary(0.0);
    covered.reservoir = Some(Reservoir {
        area_m2: AREA,
        conductance_w_k: 0.0,
        temperature_k: INITIAL_T,
    });
    let thermal = DepositionSimulation::new(
        vec![
            birth(0.25, cell(1, INITIAL_T, -1.0)),
            TimedEvent {
                time_s: 0.75,
                event: Event::Boundary(CellId(1), covered),
            },
        ],
        &[(MaterialId(1), &CALORIC)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 2,
            max_substeps: 6,
        },
    )
    .unwrap();
    if shared {
        CoatingSimulation::new_shared_chamber(
            thermal,
            &[ChamberCoatingSurface {
                cell_id: CellId(1),
                kinetics: SurfaceKinetics {
                    accommodation_coefficient: 1.0,
                    mass_transfer_coefficient_m_s: None,
                },
                initial_mass_kg: 0.0,
            }],
            &PHASE,
            &deposit_h,
            finite_bath(30.0),
            || Ok(()),
        )
        .unwrap()
    } else {
        CoatingSimulation::new(
            thermal,
            &[surface(1, 0.0, 1.0, bath(20.0), &deposit_h)],
            || Ok(()),
        )
        .unwrap()
    }
}

fn advance_history<'a>(
    initial: &CoatingSimulation<'a>,
    target: f64,
    checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<(CoatingSimulation<'a>, CoupledCoatingAudit), ThermalSublimationError> {
    if initial.shared_reservoir().is_some() {
        initial.advance_coupled_chamber_to(target, config(0.125), chamber_controls(), checkpoint)
    } else {
        initial.advance_coupled_to(target, config(0.125), phase(), checkpoint)
    }
}

#[test]
fn both_coupled_modes_preserve_birth_coverage_whole_call_budget_and_final_cancellation() {
    for shared in [false, true] {
        let initial = event_history(shared);
        let before = initial.clone();
        let mut polls = 0;
        let (next, audit) = advance_history(&initial, 1.0, || {
            polls += 1;
            Ok(())
        })
        .unwrap();
        assert_eq!(next.time_s(), 1.0);
        assert_eq!(audit.substeps, 6);
        assert_eq!(audit.phase_steps, 6);
        assert_eq!(audit.events, 2);
        assert_eq!(audit.incoming_mass_kg, SUBSTRATE_MASS);
        close(audit.external_heat_j, -0.5, 1e-14);
        assert_eq!(audit.last_steps[&CellId(1)].state.area_m2, 0.0);
        assert_eq!(audit.last_steps[&CellId(1)].phase.transferred_mass_kg, 0.0);
        assert_eq!(audit.last_steps[&CellId(1)].phase.supplied_heat_j, 0.0);
        assert_eq!(next.cells()[&CellId(1)].boundary.source_w, 0.0);
        close(audit.energy_residual_j, 0.0, 1e-6);
        assert!(audit.absolute_energy_residual_j <= 1e-6);

        let mut no_steps = initial.clone();
        no_steps.thermal.budget.max_substeps = 0;
        let (born, birth_audit) = advance_history(&no_steps, 0.25, || Ok(())).unwrap();
        assert_eq!(born.coatings()[&CellId(1)].mass_kg, 0.0);
        assert_eq!(born.shared_reservoir(), initial.shared_reservoir());
        assert_eq!(birth_audit.substeps, 0);
        assert_eq!(birth_audit.phase_steps, 0);

        let mut limited = initial.clone();
        limited.thermal.budget.max_substeps = 5;
        assert!(matches!(
            advance_history(&limited, 1.0, || Ok(())),
            Err(ThermalSublimationError::Thermal(DepositionError::Budget(_)))
        ));
        same_state(&limited, &before);
        assert_eq!(limited.shared_reservoir(), before.shared_reservoir());

        let mut visits = 0;
        assert!(matches!(
            advance_history(&initial, 1.0, || {
                visits += 1;
                if visits == polls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(ThermalSublimationError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(visits, polls);
        same_state(&initial, &before);
        assert_eq!(initial.shared_reservoir(), before.shared_reservoir());
    }
}
