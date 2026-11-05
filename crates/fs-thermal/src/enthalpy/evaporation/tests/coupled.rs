//! Heat, latent transfer and changing inventory on the original physical clock.
use super::*;

// Independent evaluation of the analytic fixture's Stefan-film mass flux.
// The production constitutive/flux routines are deliberately not the oracle.
fn analytic_flux(temperature_k: f64, gas: GasFilm) -> f64 {
    let gas_constant = 8.314_462_618_153_24;
    let molar_mass = 0.018;
    let beta = LATENT * molar_mass / gas_constant;
    let saturation = 3600.0 * (beta * (1.0 / 300.0 - 1.0 / temperature_k)).exp();
    gas.mass_transfer_coefficient_m_s * gas.total_pressure_pa * molar_mass
        / (gas_constant * gas.temperature_k)
        * ((gas.total_pressure_pa - gas.vapor_partial_pressure_pa)
            / (gas.total_pressure_pa - saturation))
            .ln()
}

fn strong_film() -> GasFilm {
    GasFilm {
        mass_transfer_coefficient_m_s: 2.0,
        ..film(0.0)
    }
}

fn narrow_phase() -> SolverConfig {
    SolverConfig {
        specific_enthalpy_bracket_j_kg: [-10_000.0, 1000.0],
        ..phase()
    }
}

fn same_timeline(actual: &DepositionSimulation<'_>, expected: &DepositionSimulation<'_>) {
    assert_eq!(actual.time_s, expected.time_s);
    assert_eq!(actual.cursor, expected.cursor);
    assert_eq!(actual.cells, expected.cells);
    assert_eq!(actual.links, expected.links);
    assert_eq!(actual.covered_areas, expected.covered_areas);
    assert_eq!(actual.schedule, expected.schedule);
    assert_eq!(actual.budget, expected.budget);
}

#[test]
fn strong_heating_feeds_latent_transfer_without_a_sensible_only_overshoot() {
    let gas = strong_film();
    let mass_rate = analytic_flux(300.0, gas) * 0.01;
    let source_w = mass_rate * LATENT;
    let thermal = initial(source_w);
    let original = thermal.cells()[&CellId(1)].clone();
    assert!(source_w / original.mass_kg > narrow_phase().specific_enthalpy_bracket_j_kg[1]);
    assert!(matches!(
        thermal.advance_evaporating_to(
            CellId(1),
            &Liquid,
            gas,
            1.0,
            config(1.0),
            narrow_phase(),
            || Ok(()),
        ),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));

    let (next, audit) = thermal
        .advance_coupled_evaporating_to(
            CellId(1),
            &Liquid,
            gas,
            1.0,
            config(1.0),
            narrow_phase(),
            || Ok(()),
        )
        .unwrap();
    close(next.temperature_k(CellId(1)).unwrap(), 300.0, 1e-10);
    close(next.mass_kg(), original.mass_kg - mass_rate, 1e-17);
    close(audit.transferred_mass_kg, mass_rate, 1e-17);
    close(audit.latent_energy_j, source_w, 1e-8);
    close(audit.vapor_energy_j, source_w, 1e-8);
    close(audit.external_heat_j, source_w, 1e-10);
    close(audit.last_step.unwrap().supplied_heat_j, source_w, 1e-10);
    close(audit.energy_residual_j, 0.0, 1e-8);
    close(audit.mass_balance_residual_kg, 0.0, 1e-17);
    assert_eq!(audit.substeps, 1);
    assert_eq!(audit.phase_steps, 1);
    assert_eq!(next.time_s(), 1.0);
    assert_eq!(next.cells()[&CellId(1)].geometry, original.geometry);
    assert_eq!(next.cells()[&CellId(1)].boundary, original.boundary);
}

#[test]
fn contact_convection_reservoir_and_source_share_one_conservative_heat_interval() {
    let mut liquid = cell(1, 0.0, 2.0);
    liquid.boundary.ambient_k = 310.0;
    liquid.boundary.htc_w_m2_k = 10.0;
    liquid.boundary.reservoir = Some(Reservoir {
        area_m2: 0.002,
        conductance_w_k: 0.2,
        temperature_k: 320.0,
    });
    let thermal = scheduled(
        vec![
            birth(0.0, liquid),
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
        10,
    )
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0;
    let duration = 0.25;
    let area = 0.007;
    let gas = film(0.0);
    // Start-state powers: 3 W contact + 0.7 W convection + 4 W reservoir + 2 W source.
    let supplied_heat = 9.7 * duration;
    let (mut lo, mut hi) = (-20.0 * CP, 30.0 * CP);
    for _ in 0..80 {
        let h = lo + 0.5 * (hi - lo);
        let dm = analytic_flux(300.0 + h / CP, gas) * area * duration;
        if 0.02 * h + dm * LATENT > supplied_heat {
            hi = h;
        } else {
            lo = h;
        }
    }
    let expected_h = lo + 0.5 * (hi - lo);
    let expected_transfer = analytic_flux(300.0 + expected_h / CP, gas) * area * duration;

    let (next, audit) = thermal
        .advance_coupled_evaporating_to(
            CellId(1),
            &Liquid,
            gas,
            duration,
            config(duration),
            phase(),
            || Ok(()),
        )
        .unwrap();
    close(next.cells()[&CellId(1)].specific_h_j_kg, expected_h, 1e-6);
    close(audit.transferred_mass_kg, expected_transfer, 1e-14);
    close(
        next.cells()[&CellId(2)].specific_h_j_kg,
        10.0 * CP - 0.75 / 0.02,
        1e-10,
    );
    close(
        audit.last_step.as_ref().unwrap().supplied_heat_j,
        2.425,
        1e-13,
    );
    close(
        audit.last_step.as_ref().unwrap().liquid.area_m2,
        area,
        1e-16,
    );
    close(audit.external_heat_j, 1.675, 1e-13);
    close(
        audit.vapor_energy_j,
        expected_transfer * (expected_h + LATENT),
        1e-8,
    );
    close(
        next.energy_j() - thermal.energy_j() + audit.vapor_energy_j,
        1.675,
        1e-8,
    );
    close(audit.energy_residual_j, 0.0, 1e-8);
    close(
        next.mass_kg() + audit.transferred_mass_kg,
        thermal.mass_kg(),
        1e-16,
    );
    assert_eq!(audit.substeps, 1);
    assert_eq!(audit.phase_steps, 1);
    assert_eq!(next.links(), thermal.links());
}

#[test]
fn evaporated_mass_tightens_the_original_cfl_and_shares_the_step_budget() {
    let gas = strong_film();
    let rate = analytic_flux(300.0, gas) * 0.009;
    let mut liquid = cell(1, 0.0, rate * LATENT);
    liquid.boundary.reservoir = Some(Reservoir {
        area_m2: 0.001,
        conductance_w_k: 100.0,
        temperature_k: 300.0,
    });
    let thermal = scheduled(vec![birth(0.0, liquid)], 2)
        .advance_to(0.0, config(1.0), || Ok(()))
        .unwrap()
        .0;
    let before = thermal.clone();
    let first_bound = 0.8 / (100.0 * (1.0 / CP) / 0.02);
    let target = 2.0 * first_bound;
    let next_mass = 0.02 - rate * first_bound;
    let second_bound = 0.8 / (100.0 * (1.0 / CP) / next_mass);
    assert!(second_bound < first_bound);
    // Two steps would suffice if the original mass were incorrectly retained.
    assert!(matches!(
        thermal.advance_coupled_evaporating_to(
            CellId(1),
            &Liquid,
            gas,
            target,
            config(target),
            narrow_phase(),
            || Ok(()),
        ),
        Err(ThermalEvaporationError::Thermal(DepositionError::Budget(_)))
    ));
    same_timeline(&thermal, &before);

    let mut admitted = thermal.clone();
    admitted.budget.max_substeps = 3;
    let (next, audit) = admitted
        .advance_coupled_evaporating_to(
            CellId(1),
            &Liquid,
            gas,
            target,
            config(target),
            narrow_phase(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(audit.substeps, 3);
    assert_eq!(audit.phase_steps, 3);
    close(next.temperature_k(CellId(1)).unwrap(), 300.0, 1e-10);
    close(next.mass_kg(), 0.02 - rate * target, 1e-16);
    close(audit.external_heat_j, rate * LATENT * target, 1e-8);
    close(audit.energy_residual_j, 0.0, 1e-8);
    assert_eq!(next.time_s(), target);
}

#[test]
fn birth_and_cover_events_clip_both_surface_transfer_and_authored_heat() {
    let gas = strong_film();
    let rate = analytic_flux(300.0, gas) * 0.01;
    let source = rate * LATENT;
    let mut covered = boundary(0.0);
    covered.reservoir = Some(Reservoir {
        area_m2: 0.01,
        conductance_w_k: 0.0,
        temperature_k: 300.0,
    });
    let thermal = scheduled(
        vec![
            birth(0.5, cell(1, 0.0, source)),
            TimedEvent {
                time_s: 1.0,
                event: Event::Boundary(CellId(1), covered),
            },
        ],
        4,
    );
    let (born, birth_audit) = thermal
        .advance_coupled_evaporating_to(
            CellId(1),
            &Liquid,
            gas,
            0.5,
            config(1.0),
            narrow_phase(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(born.mass_kg(), 0.02);
    assert_eq!(birth_audit.phase_steps, 0);
    assert_eq!(birth_audit.substeps, 0);
    assert_eq!(birth_audit.external_heat_j, 0.0);

    let (next, audit) = thermal
        .advance_coupled_evaporating_to(
            CellId(1),
            &Liquid,
            gas,
            2.0,
            config(1.0),
            narrow_phase(),
            || Ok(()),
        )
        .unwrap();
    close(next.mass_kg(), 0.02 - 0.5 * rate, 1e-17);
    close(next.temperature_k(CellId(1)).unwrap(), 300.0, 1e-10);
    close(audit.external_heat_j, 0.5 * source, 1e-8);
    close(audit.vapor_energy_j, 0.5 * source, 1e-8);
    close(audit.energy_residual_j, 0.0, 1e-8);
    assert_eq!(audit.incoming_mass_kg, 0.02);
    assert_eq!(audit.incoming_energy_j, 0.0);
    assert_eq!(audit.events, 2);
    assert_eq!(audit.substeps, 2);
    assert_eq!(audit.phase_steps, 2);
    assert_eq!(next.exposed_area_m2(CellId(1)), Ok(0.0));
    assert_eq!(next.cells()[&CellId(1)].boundary, covered);
    assert_eq!(audit.last_step.as_ref().unwrap().transferred_mass_kg, 0.0);
    assert_eq!(audit.last_step.as_ref().unwrap().supplied_heat_j, 0.0);
    assert_eq!(next.time_s(), 2.0);
}

#[test]
fn relative_balance_survives_shifted_birth_enthalpy_and_large_unchanged_neighbor_energy() {
    struct ShiftedCaloric(f64);
    impl ThermalLaw for ShiftedCaloric {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            CALORIC.temperature_k(h - self.0)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0 / CP
        }
    }
    let gas = film(0.0);
    let duration = 0.5;
    // This exactly representable increment makes chart quantization identical
    // on every datum while exercising both sensible storage and latent export.
    let increment = 512.0;
    let transfer = analytic_flux(300.0 + increment / CP, gas) * 0.01 * duration;
    let heat = 0.02 * increment + transfer * LATENT;
    let large_datum = 2.0_f64.powi(40);
    let neighbor_datum = 2.0_f64.powi(80);
    let mut reference = None;
    for datum in [0.0, large_datum, -large_datum] {
        let liquid_law = ShiftedCaloric(datum);
        let neighbor_law = ShiftedCaloric(neighbor_datum);
        let mut neighbor = cell(2, neighbor_datum, 0.0);
        neighbor.material_id = MaterialId(2);
        let thermal = DepositionSimulation::new(
            vec![
                birth(0.0, neighbor.clone()),
                birth(0.5, cell(1, datum, heat / duration)),
            ],
            &[
                (MaterialId(1), &liquid_law as &dyn ThermalLaw),
                (MaterialId(2), &neighbor_law as &dyn ThermalLaw),
            ],
            DepositionBudget {
                max_cells: 2,
                max_links: 0,
                max_events: 2,
                max_substeps: 2,
            },
        )
        .unwrap()
        .advance_to(0.0, config(1.0), || Ok(()))
        .unwrap()
        .0;
        let (next, audit) = thermal
            .advance_coupled_evaporating_to(
                CellId(1),
                &Liquid,
                gas,
                1.0,
                config(1.0),
                SolverConfig {
                    specific_enthalpy_bracket_j_kg: [datum - 1024.0, datum + 1024.0],
                    ..phase()
                },
                || Ok(()),
            )
            .unwrap();
        assert_eq!(next.cells()[&CellId(2)], neighbor);
        assert_eq!(next.cells()[&CellId(1)].specific_h_j_kg - datum, increment);
        assert!(audit.transferred_mass_kg > 0.0);
        close(audit.transferred_mass_kg, transfer, 1e-17);
        close(audit.external_heat_j, heat, 1e-12);
        close(audit.energy_residual_j, 0.0, 1e-8);
        assert_eq!(audit.incoming_mass_kg, 0.02);
        assert_eq!(audit.incoming_energy_j, 0.02 * datum);
        assert_eq!(audit.events, 1);
        assert_eq!(audit.substeps, 2);
        assert_eq!(audit.phase_steps, 1);
        let physical = (
            next.cells()[&CellId(1)].mass_kg,
            audit.transferred_mass_kg,
            audit.latent_energy_j,
            audit.external_heat_j,
            audit.mass_balance_residual_kg,
            audit.energy_residual_j,
        );
        if let Some(expected) = reference {
            assert_eq!(physical, expected);
        } else {
            reference = Some(physical);
        }
        if datum == 0.0 {
            // Aggregate energy subtraction loses the entire sensible change.
            // It cannot replace the admitted relative open-system balance.
            assert_eq!(next.energy_j(), thermal.energy_j());
            let rounded_balance = next.energy_j() - thermal.energy_j() + audit.vapor_energy_j
                - audit.incoming_energy_j
                - audit.external_heat_j;
            assert!(rounded_balance.abs() > 1.0);
            assert!(audit.stored_energy_change_j > 10.0);
        }
    }
}

#[test]
fn later_event_failure_and_final_cancellation_leave_the_timeline_unchanged() {
    let mut invalid = boundary(0.0);
    invalid.reservoir = Some(Reservoir {
        area_m2: 0.02,
        conductance_w_k: 0.0,
        temperature_k: 300.0,
    });
    let thermal = scheduled(
        vec![
            birth(0.0, cell(1, 0.0, 2.0)),
            TimedEvent {
                time_s: 0.5,
                event: Event::Boundary(CellId(1), invalid),
            },
        ],
        10,
    );
    let before = thermal.clone();
    let (prefix, prefix_audit) = thermal
        .advance_coupled_evaporating_to(
            CellId(1),
            &Liquid,
            film(0.0),
            0.25,
            config(0.25),
            phase(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(prefix_audit.phase_steps, 1);
    assert!(prefix.mass_kg() < 0.02);
    assert!(matches!(
        thermal.advance_coupled_evaporating_to(
            CellId(1),
            &Liquid,
            film(0.0),
            1.0,
            config(0.25),
            phase(),
            || Ok(()),
        ),
        Err(ThermalEvaporationError::Thermal(
            DepositionError::InvalidInput(_)
        ))
    ));
    same_timeline(&thermal, &before);

    let thermal = initial(2.0);
    let before = thermal.clone();
    let mut checkpoints = 0;
    let (_, accepted) = thermal
        .advance_coupled_evaporating_to(
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
    assert_eq!(accepted.phase_steps, 4);
    let mut visited = 0;
    let cancelled = thermal.advance_coupled_evaporating_to(
        CellId(1),
        &Liquid,
        film(0.0),
        1.0,
        config(0.25),
        phase(),
        || {
            visited += 1;
            if visited == checkpoints {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        cancelled,
        Err(ThermalEvaporationError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(visited, checkpoints);
    same_timeline(&thermal, &before);
}
