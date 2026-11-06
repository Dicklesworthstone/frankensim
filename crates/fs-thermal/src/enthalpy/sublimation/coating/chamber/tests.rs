//! G0/G3: one finite gas supply across native thermal events and supported cells.

use super::super::ChamberCoatingSurface;
use super::*;
use crate::sublimation::chamber::ChamberSolverConfig;

fn chamber_config() -> ChamberSolverConfig {
    ChamberSolverConfig {
        surface_solver: phase(),
        relative_mass_tolerance: 1e-11,
        absolute_mass_tolerance_kg: 1e-16,
        max_iterations: 128,
        max_surfaces: 4,
    }
}

fn chamber_bath(pressure_pa: f64) -> VaporReservoir {
    VaporReservoir::Finite {
        vapor_mass_kg: pressure_pa * 0.001 * 0.018 / (MOLAR_GAS_CONSTANT * INITIAL_T),
        volume_m3: 0.001,
        temperature_k: INITIAL_T,
    }
}

fn gas_mass(reservoir: VaporReservoir) -> f64 {
    match reservoir {
        VaporReservoir::Finite { vapor_mass_kg, .. } => vapor_mass_kg,
        _ => panic!("shared gas must remain finite"),
    }
}

fn selected(id: u64, initial_mass_kg: f64) -> ChamberCoatingSurface {
    ChamberCoatingSurface {
        cell_id: CellId(id),
        kinetics: SurfaceKinetics {
            accommodation_coefficient: 1.0,
            mass_transfer_coefficient_m_s: None,
        },
        initial_mass_kg,
    }
}

fn shared<'a>(
    events: Vec<TimedEvent>,
    surfaces: &[ChamberCoatingSurface],
    deposit_enthalpy: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
) -> CoatingSimulation<'a> {
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(1), &CALORIC)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 16,
            max_substeps: 64,
        },
    )
    .unwrap();
    CoatingSimulation::new_shared_chamber(
        thermal,
        surfaces,
        &PHASE,
        deposit_enthalpy,
        chamber_bath(30.0),
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn two_bare_substrates_compete_for_one_inventory_and_replay_in_cell_order() {
    let events = vec![
        birth(0.0, cell(1, INITIAL_T, 0.0)),
        birth(0.0, cell(2, INITIAL_T, 0.0)),
    ];
    let initial = shared(
        events.clone(),
        &[selected(1, 0.0), selected(2, 0.0)],
        &deposit_h,
    );
    let reordered = shared(events, &[selected(2, 0.0), selected(1, 0.0)], &deposit_h);
    let (next, audit) = initial
        .advance_chamber_to(0.5, config(0.5), chamber_config(), || Ok(()))
        .unwrap();
    let (other, other_audit) = reordered
        .advance_chamber_to(0.5, config(0.5), chamber_config(), || Ok(()))
        .unwrap();
    let original_vapor = gas_mass(initial.shared_reservoir().unwrap());
    let vapor = gas_mass(next.shared_reservoir().unwrap());
    assert!(next.deposit_mass_kg() > 0.0);
    assert!(next.deposit_mass_kg() < original_vapor);
    close(next.deposit_mass_kg() + vapor, original_vapor, 1e-18);
    assert_eq!(
        next.coatings()[&CellId(1)].mass_kg,
        next.coatings()[&CellId(2)].mass_kg
    );
    for coating in next.coatings().values() {
        assert_eq!(Some(coating.reservoir), next.shared_reservoir());
    }
    assert!(audit.pressure_iterations > 0);
    close(
        audit.stored_energy_change_j - audit.thermal.incoming_energy_j + audit.vapor_energy_j,
        0.0,
        1e-6,
    );
    closes(&audit);
    same_state(&next, &other);
    assert_eq!(next.shared_reservoir(), other.shared_reservoir());
    assert_eq!(audit, other_audit);
}

#[test]
fn late_substrate_joins_current_vapor_without_resetting_the_chamber() {
    let initial = shared(
        vec![
            birth(0.0, cell(1, INITIAL_T, 0.0)),
            birth(0.5, cell(2, INITIAL_T, 0.0)),
        ],
        &[selected(1, 0.0), selected(2, 0.0)],
        &deposit_h,
    );
    let (middle, first) = initial
        .advance_chamber_to(0.5, config(0.25), chamber_config(), || Ok(()))
        .unwrap();
    assert_eq!(middle.coatings()[&CellId(2)].mass_kg, 0.0);
    assert_eq!(
        Some(middle.coatings()[&CellId(2)].reservoir),
        middle.shared_reservoir()
    );
    assert!(
        gas_mass(middle.shared_reservoir().unwrap())
            < gas_mass(initial.shared_reservoir().unwrap())
    );
    let (resumed, second) = middle
        .advance_chamber_to(1.0, config(0.25), chamber_config(), || Ok(()))
        .unwrap();
    let (whole, all) = initial
        .advance_chamber_to(1.0, config(0.25), chamber_config(), || Ok(()))
        .unwrap();
    assert!(resumed.coatings()[&CellId(2)].mass_kg > 0.0);
    close(
        resumed.deposit_mass_kg() + gas_mass(resumed.shared_reservoir().unwrap()),
        gas_mass(initial.shared_reservoir().unwrap()),
        1e-18,
    );
    same_state(&resumed, &whole);
    assert_eq!(resumed.shared_reservoir(), whole.shared_reservoir());
    assert_eq!(
        first.thermal.events + second.thermal.events,
        all.thermal.events
    );
    closes(&first);
    closes(&second);
    closes(&all);
}

#[test]
fn no_active_substrates_preserves_gas_until_birth_and_endpoint_coverage_blocks_exchange() {
    let link = ThermalLink {
        id: 9,
        a: CellId(1),
        b: CellId(2),
        area_m2: AREA,
        conductance_w_k: 0.0,
    };
    let initial = shared(
        vec![
            birth(0.5, cell(1, INITIAL_T, 0.0)),
            birth(0.5, cell(2, INITIAL_T, 0.0)),
            TimedEvent {
                time_s: 0.5,
                event: Event::Contact(link),
            },
        ],
        &[selected(1, 0.0), selected(2, 0.0)],
        &deposit_h,
    );
    let (waiting, before) = initial
        .advance_chamber_to(0.25, config(0.25), chamber_config(), || Ok(()))
        .unwrap();
    assert!(waiting.coatings().is_empty());
    assert_eq!(waiting.shared_reservoir(), initial.shared_reservoir());
    assert_eq!(before.phase_steps, 0);
    let (covered, after) = waiting
        .advance_chamber_to(1.0, config(0.25), chamber_config(), || Ok(()))
        .unwrap();
    assert_eq!(covered.deposit_mass_kg(), 0.0);
    assert_eq!(covered.shared_reservoir(), initial.shared_reservoir());
    assert_eq!(after.transferred_mass_kg, 0.0);
    closes(&after);
}

#[test]
fn chamber_heat_uses_both_material_capacities_and_large_substrate_datum_is_independent() {
    let mut inert = selected(1, 0.005);
    inert.kinetics.accommodation_coefficient = 0.0;
    let events = vec![birth(0.0, cell(1, INITIAL_T, 10.0))];
    let initial = shared(events, &[inert], &deposit_h);
    let (heated, audit) = initial
        .advance_chamber_to(0.25, config(0.25), chamber_config(), || Ok(()))
        .unwrap();
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        INITIAL_T + 2.5 / 20.0,
        1e-9,
    );
    close(audit.thermal.external_heat_j, 2.5, 1e-14);
    close(
        audit.stored_energy_change_j - audit.thermal.incoming_energy_j,
        2.5,
        1e-8,
    );
    assert_eq!(heated.shared_reservoir(), initial.shared_reservoir());
    closes(&audit);

    let shifted = Caloric { datum: 1e9 };
    let mut substrate = cell(1, INITIAL_T, 10.0);
    substrate.specific_h_j_kg += shifted.datum;
    let thermal = DepositionSimulation::new(
        vec![birth(0.0, substrate)],
        &[(MaterialId(1), &shifted)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 2,
            max_substeps: 8,
        },
    )
    .unwrap();
    let shifted_state = CoatingSimulation::new_shared_chamber(
        thermal,
        &[inert],
        &PHASE,
        &deposit_h,
        chamber_bath(30.0),
        || Ok(()),
    )
    .unwrap();
    let (shifted_heat, shifted_audit) = shifted_state
        .advance_chamber_to(0.25, config(0.25), chamber_config(), || Ok(()))
        .unwrap();
    close(
        shifted_heat.temperature_k(CellId(1)).unwrap(),
        heated.temperature_k(CellId(1)).unwrap(),
        1e-8,
    );
    closes(&shifted_audit);
}

#[test]
fn invalid_shared_controls_late_cancellation_and_individual_replacement_preserve_state() {
    let initial = shared(
        vec![
            birth(0.0, cell(1, INITIAL_T, 0.0)),
            birth(0.5, cell(2, INITIAL_T, 0.0)),
        ],
        &[selected(1, 0.0), selected(2, 0.0)],
        &deposit_h,
    );
    let (mut ready, _) = initial
        .advance_chamber_to(0.0, config(0.25), chamber_config(), || Ok(()))
        .unwrap();
    let original = ready.clone();
    assert!(
        ready
            .set_reservoir(CellId(1), chamber_bath(10.0), || Ok(()))
            .is_err()
    );
    assert!(
        ready
            .advance_to(0.5, config(0.25), phase(), || Ok(()))
            .is_err()
    );
    let mut invalid = chamber_config();
    invalid.max_iterations = 0;
    assert!(
        ready
            .advance_chamber_to(0.0, config(0.25), invalid, || Ok(()))
            .is_err()
    );
    invalid = chamber_config();
    invalid.max_surfaces = 1;
    assert!(
        ready
            .advance_chamber_to(0.5, config(0.25), invalid, || Ok(()))
            .is_err()
    );

    let mut calls = 0usize;
    ready
        .advance_chamber_to(1.0, config(0.25), chamber_config(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let mut attempted = 0usize;
    assert!(
        ready
            .advance_chamber_to(1.0, config(0.25), chamber_config(), || {
                attempted += 1;
                if attempted == calls - 1 {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    same_state(&ready, &original);
    assert_eq!(ready.shared_reservoir(), original.shared_reservoir());

    let replacement = chamber_bath(15.0);
    ready.set_shared_reservoir(replacement, || Ok(())).unwrap();
    assert_eq!(ready.shared_reservoir(), Some(replacement));
    assert_eq!(ready.coatings()[&CellId(1)].reservoir, replacement);
    let (born, _) = ready
        .advance_chamber_to(0.5, config(0.25), chamber_config(), || Ok(()))
        .unwrap();
    assert_eq!(
        Some(born.coatings()[&CellId(2)].reservoir),
        born.shared_reservoir()
    );
}

fn ready_for_supplied_chamber_heat() -> CoatingSimulation<'static> {
    shared(
        vec![
            birth(0.0, cell(1, INITIAL_T, 11.0)),
            birth(0.0, cell(2, INITIAL_T, -3.0)),
            birth(1.0, cell(3, INITIAL_T, 0.0)),
        ],
        &[selected(1, 1e-3), selected(2, 0.0), selected(3, 0.0)],
        &deposit_h,
    )
    .advance_chamber_to(0.0, config(0.1), chamber_config(), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn supplied_chamber_heat_uses_ids_and_keeps_unheated_bare_receivers() {
    let original = ready_for_supplied_chamber_heat();
    let mut next = original.clone();
    let (ids, step) = next
        .exchange_chamber_with_heat(0.1, &[(CellId(1), 120.0)], chamber_config(), || Ok(()))
        .unwrap();
    assert_eq!(ids, vec![CellId(1), CellId(2)]);
    assert_eq!(step.supplied_heat_j, 120.0);
    assert_eq!(step.surfaces[0].phase.supplied_heat_j, 120.0);
    assert_eq!(step.surfaces[1].phase.supplied_heat_j, 0.0);
    assert!(next.coatings()[&CellId(2)].mass_kg > 0.0);
    assert!(next.temperature_k(CellId(1)).unwrap() > next.temperature_k(CellId(2)).unwrap());
    for id in &ids {
        assert_eq!(next.cells()[id].mass_kg, original.cells()[id].mass_kg);
        assert_eq!(next.cells()[id].boundary, original.cells()[id].boundary);
        assert_eq!(Some(next.coatings()[id].reservoir), next.shared_reservoir());
    }
    assert_eq!(next.time_s(), original.time_s());
    assert_eq!(next.pending_events(), original.pending_events());
    close(
        next.mass_kg() + gas_mass(next.shared_reservoir().unwrap()),
        original.mass_kg() + gas_mass(original.shared_reservoir().unwrap()),
        1e-16,
    );
    close(
        next.energy_j() - original.energy_j() + step.vapor_energy_j - step.supplied_heat_j,
        step.energy_residual_j,
        1e-7,
    );
    close(step.energy_residual_j, 0.0, 1e-7);
    let mut explicit_zero = original;
    let (other_ids, other_step) = explicit_zero
        .exchange_chamber_with_heat(
            0.1,
            &[(CellId(2), 0.0), (CellId(1), 120.0)],
            chamber_config(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(other_ids, ids);
    assert_eq!(other_step, step);
    same_state(&next, &explicit_zero);
    assert_eq!(next.shared_reservoir(), explicit_zero.shared_reservoir());
}

#[test]
fn supplied_shared_heat_refusal_and_final_cancellation_preserve_every_owner() {
    let original = ready_for_supplied_chamber_heat();
    for invalid in [
        &[(CellId(1), 1.0), (CellId(1), 2.0)][..],
        &[(CellId(99), 1.0)][..],
        &[(CellId(3), 1.0)][..],
        &[(CellId(1), f64::NAN)][..],
        &[(CellId(2), f64::INFINITY)][..],
    ] {
        let mut next = original.clone();
        assert!(
            next.exchange_chamber_with_heat(0.1, invalid, chamber_config(), || Ok(()))
                .is_err()
        );
        same_state(&next, &original);
        assert_eq!(next.shared_reservoir(), original.shared_reservoir());
    }
    let mut restricted = chamber_config();
    restricted.max_surfaces = 1;
    let mut next = original.clone();
    assert!(
        next.exchange_chamber_with_heat(0.1, &[(CellId(1), 120.0)], restricted, || Ok(()))
            .is_err()
    );
    same_state(&next, &original);
    assert_eq!(next.shared_reservoir(), original.shared_reservoir());

    let mut probe = original.clone();
    let mut checkpoints = 0;
    probe
        .exchange_chamber_with_heat(0.1, &[(CellId(1), 120.0)], chamber_config(), || {
            checkpoints += 1;
            Ok(())
        })
        .unwrap();
    let mut calls = 0;
    let result =
        next.exchange_chamber_with_heat(0.1, &[(CellId(1), 120.0)], chamber_config(), || {
            calls += 1;
            if calls == checkpoints {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
    assert_eq!(
        result,
        Err(ThermalSublimationError::Thermal(DepositionError::Cancelled))
    );
    same_state(&next, &original);
    assert_eq!(next.shared_reservoir(), original.shared_reservoir());
}
