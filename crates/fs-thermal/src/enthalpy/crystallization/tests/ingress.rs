//! Original incoming crystal inventory, energy and event-time semantics.

use super::*;
use crate::crystallization::CrystallizationAdaptivity;
use crate::enthalpy::crystallization::implicit::CrystallizationImplicitControls;

fn cell(id: u64, h: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(3),
        mass_kg: 0.01,
        specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64 + 0.01, 0.0, 0.0],
            volume_m3: 1e-5,
            surface_area_m2: 0.01,
        },
        boundary: boundary(0.0),
    }
}
fn birth(time_s: f64, cell: DepositedCell) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Birth(cell),
    }
}
fn feed(events: Vec<TimedEvent>) -> DepositionSimulation<'static> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(3), &AMORPHOUS), (MaterialId(4), &MOLD)],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 32,
            max_substeps: 10_000,
        },
    )
    .unwrap()
}
fn crystal_mass(simulation: &CrystallizationSimulation<'_>) -> f64 {
    simulation
        .states()
        .iter()
        .map(|(id, state)| simulation.cells()[id].mass_kg * state.crystalline_mass_fraction())
        .sum()
}

#[test]
fn active_history_keeps_its_original_heat_and_crystal_inventory() {
    let law = law();
    let thermal = finite_bath();
    let h = thermal.cells()[&CellId(10)].specific_h_j_kg;
    let incoming = law
        .state_at_enthalpy(h, 0.6, controls(0.125).caloric, || false)
        .unwrap();
    let initial = CrystallizationSimulation::new_with_incoming_states(
        thermal.clone(),
        &[material(&law)],
        &[(CellId(10), incoming)],
        controls(0.125).caloric,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(initial.states()[&CellId(10)], incoming);
    assert_eq!(initial.energy_j(), thermal.energy_j());
    assert_eq!(initial.mass_kg(), thermal.mass_kg());
    assert_eq!(
        initial.temperature_k(CellId(10)).unwrap(),
        incoming.temperature_k()
    );
    assert!(
        initial.temperature_k(CellId(10)).unwrap() > thermal.temperature_k(CellId(10)).unwrap()
    );
    assert_eq!(
        initial.temperature_k(CellId(11)).unwrap(),
        thermal.temperature_k(CellId(11)).unwrap()
    );
    let (_, no_elapsed_time) = initial.advance_to(0.0, controls(0.125), || Ok(())).unwrap();
    assert_eq!(no_elapsed_time.crystallized_mass_kg, 0.0);
    assert_eq!(no_elapsed_time.thermal.incoming_mass_kg, 0.0);

    let (next, audit) = initial
        .advance_to(0.25, controls(0.125), || Ok(()))
        .unwrap();
    assert!(next.states()[&CellId(10)].exposure() > incoming.exposure());
    close(
        audit.crystallized_mass_kg,
        crystal_mass(&next) - crystal_mass(&initial),
        1e-16,
    );
    assert_eq!(audit.thermal.incoming_mass_kg, 0.0);
    assert_eq!(audit.thermal.incoming_energy_j, 0.0);
    close(next.energy_j(), initial.energy_j(), 1e-8);
    assert!(audit.absolute_energy_residual_j < 1e-7);
    let (replayed, replay_audit) = initial
        .clone()
        .advance_to(0.25, controls(0.125), || Ok(()))
        .unwrap();
    assert_eq!(replayed.states(), next.states());
    assert_eq!(replay_audit, audit);
}

#[test]
fn cold_precrystallized_births_keep_energy_and_do_not_age_before_activation() {
    let law = law();
    let cold = law.state_at_temperature(280.0, 0.9).unwrap();
    let later = law.state_at_temperature(390.0, 0.3).unwrap();
    // This h is valid for the supplied mixture but outside the parent chart.
    assert!(ThermalLaw::temperature_k(&AMORPHOUS, cold.specific_enthalpy_j_kg()).is_err());
    let thermal = feed(vec![
        birth(0.0, cell(10, 50_000.0)),
        birth(0.5, cell(11, cold.specific_enthalpy_j_kg())),
        birth(1.0, cell(12, later.specific_enthalpy_j_kg())),
    ]);
    let initial = CrystallizationSimulation::new_with_incoming_states(
        thermal.clone(),
        &[material(&law)],
        &[(CellId(11), cold), (CellId(12), later)],
        controls(0.125).caloric,
        || Ok(()),
    )
    .unwrap();
    let (start, _) = initial.advance_to(0.0, controls(0.125), || Ok(())).unwrap();
    assert_eq!(start.states()[&CellId(10)].exposure(), 0.0);
    let (middle, first_audit) = initial.advance_to(0.5, controls(0.125), || Ok(())).unwrap();
    assert!(middle.states()[&CellId(10)].exposure() > 0.0);
    assert_eq!(middle.states()[&CellId(11)], cold);
    assert!(!middle.states().contains_key(&CellId(12)));
    close(
        first_audit.crystallized_mass_kg,
        crystal_mass(&middle) - 0.01 * cold.crystalline_mass_fraction(),
        1e-16,
    );
    close(first_audit.thermal.incoming_mass_kg, 0.02, 1e-16);
    close(
        first_audit.thermal.incoming_energy_j,
        0.01 * (50_000.0 + cold.specific_enthalpy_j_kg()),
        1e-10,
    );
    close(
        middle.energy_j(),
        first_audit.thermal.incoming_energy_j,
        1e-8,
    );

    let (next, second_audit) = middle.advance_to(1.0, controls(0.125), || Ok(())).unwrap();
    assert_eq!(next.states()[&CellId(12)], later);
    assert!(next.states()[&CellId(11)].exposure() > cold.exposure());
    assert!(next.states()[&CellId(10)].exposure() > middle.states()[&CellId(10)].exposure());
    close(
        second_audit.crystallized_mass_kg,
        crystal_mass(&next) - crystal_mass(&middle) - 0.01 * later.crystalline_mass_fraction(),
        1e-16,
    );
    close(second_audit.thermal.incoming_mass_kg, 0.01, 1e-16);
    close(
        second_audit.thermal.incoming_energy_j,
        0.01 * later.specific_enthalpy_j_kg(),
        1e-10,
    );
    let (whole, whole_audit) = initial.advance_to(1.0, controls(0.125), || Ok(())).unwrap();
    assert_eq!(whole.states(), next.states());
    close(
        whole_audit.crystallized_mass_kg,
        first_audit.crystallized_mass_kg + second_audit.crystallized_mass_kg,
        1e-16,
    );
    // Original zero-exposure behavior still refuses the same unsupported feed.
    let default = CrystallizationSimulation::new(thermal, &[material(&law)], || Ok(())).unwrap();
    assert!(matches!(
        default.advance_to(0.5, controls(0.125), || Ok(())),
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
}

#[test]
fn implicit_candidates_retain_future_crystal_feed() {
    let law = law();
    let first = law.state_at_temperature(280.0, 0.9).unwrap();
    let later = law.state_at_temperature(290.0, 0.8).unwrap();
    let initial = CrystallizationSimulation::new_with_incoming_states(
        feed(vec![
            birth(0.0, cell(10, first.specific_enthalpy_j_kg())),
            birth(0.25, cell(11, later.specific_enthalpy_j_kg())),
        ]),
        &[material(&law)],
        &[(CellId(10), first), (CellId(11), later)],
        controls(0.125).caloric,
        || Ok(()),
    )
    .unwrap()
    .with_adaptive_phase_integration(MaterialId(3), CrystallizationAdaptivity::default())
    .unwrap();
    let settings = CrystallizationImplicitControls {
        phase: controls(0.125),
        ..CrystallizationImplicitControls::default()
    };
    let (next, report) = initial
        .advance_implicit_to(0.25, settings, || Ok(()))
        .unwrap();
    assert_eq!(next.states()[&CellId(11)], later);
    assert!(next.states()[&CellId(10)].exposure() > first.exposure());
    close(
        report.phase.crystallized_mass_kg,
        crystal_mass(&next)
            - 0.01 * (first.crystalline_mass_fraction() + later.crystalline_mass_fraction()),
        1e-16,
    );
    close(
        next.energy_j(),
        0.01 * (first.specific_enthalpy_j_kg() + later.specific_enthalpy_j_kg()),
        1e-8,
    );
    let (resumed, _) = next
        .advance_implicit_to(0.375, settings, || Ok(()))
        .unwrap();
    assert!(resumed.states()[&CellId(11)].exposure() > later.exposure());
}

#[test]
fn ingress_refuses_wrong_history_ambiguous_births_and_unused_keys_atomically() {
    let law = law();
    let original = finite_bath();
    let state = law
        .state_at_enthalpy(80_000.0, 0.6, controls(0.125).caloric, || false)
        .unwrap();
    let attach = |thermal, supplied: &[(CellId, CrystallizationState)]| {
        CrystallizationSimulation::new_with_incoming_states(
            thermal,
            &[material(&law)],
            supplied,
            controls(0.125).caloric,
            || Ok(()),
        )
    };
    assert!(
        attach(
            original.clone(),
            &[(CellId(10), state), (CellId(10), state)]
        )
        .is_err()
    );
    assert!(attach(original.clone(), &[(CellId(999), state)]).is_err());
    assert!(attach(original.clone(), &[(CellId(11), state)]).is_err());
    let wrong_h = law.state_at_temperature(350.0, 0.6).unwrap();
    assert!(matches!(
        attach(original.clone(), &[(CellId(10), wrong_h)]),
        Err(CrystallizationTimelineError::Material(
            CrystallizationError::StateMismatch
        ))
    ));
    let other_law = CrystallizationLaw::new(
        &AMORPHOUS,
        &CRYSTALLINE,
        &RATE,
        CrystallizationConfig {
            identity: [38; 32],
            ..law.config()
        },
    )
    .unwrap();
    let other_state = other_law
        .state_at_enthalpy(80_000.0, 0.6, controls(0.125).caloric, || false)
        .unwrap();
    assert!(matches!(
        attach(original.clone(), &[(CellId(10), other_state)]),
        Err(CrystallizationTimelineError::Material(
            CrystallizationError::StateMismatch
        ))
    ));
    let (advanced, _) = original
        .advance_to(
            0.25,
            StepConfig {
                max_dt_s: 0.125,
                energy_tolerance_j: 1e-7,
            },
            || Ok(()),
        )
        .unwrap();
    assert!(matches!(
        attach(advanced, &[(CellId(10), state)]),
        Err(CrystallizationTimelineError::Material(
            CrystallizationError::StateMismatch
        ))
    ));
    let repeated = feed(vec![
        birth(0.0, cell(10, state.specific_enthalpy_j_kg())),
        birth(0.5, cell(10, state.specific_enthalpy_j_kg())),
    ]);
    assert!(attach(repeated, &[(CellId(10), state)]).is_err());
    let over_budget = vec![(CellId(10), state); original.budget.max_cells + 1];
    assert!(matches!(
        attach(original.clone(), &over_budget),
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::Budget(_)
        ))
    ));
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(original.cells()[&CellId(10)].specific_h_j_kg, 80_000.0);
}

#[test]
fn late_original_birth_failure_and_cancellation_publish_no_feed_or_history() {
    let law = law();
    let incoming = law.state_at_temperature(280.0, 0.9).unwrap();
    let mut invalid_cell = cell(11, incoming.specific_enthalpy_j_kg());
    invalid_cell.mass_kg = 0.0;
    let initial = CrystallizationSimulation::new_with_incoming_states(
        feed(vec![
            birth(0.0, cell(10, incoming.specific_enthalpy_j_kg())),
            birth(0.25, invalid_cell),
        ]),
        &[material(&law)],
        &[(CellId(10), incoming), (CellId(11), incoming)],
        controls(0.125).caloric,
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        initial.advance_to(0.25, controls(0.125), || Ok(())),
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::InvalidInput(_)
        ))
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert!(initial.cells().is_empty());
    assert!(initial.states().is_empty());
    assert_eq!(initial.pending_events().len(), 2);
    let mut checks = 0;
    assert!(matches!(
        initial.advance_to(0.125, controls(0.125), || {
            checks += 1;
            if checks == 12 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    assert!(checks >= 12);
    assert!(initial.states().is_empty());
    let (activated, _) = initial.advance_to(0.0, controls(0.125), || Ok(())).unwrap();
    assert_eq!(activated.states()[&CellId(10)], incoming);
    assert_eq!(activated.pending_events().len(), 1);
}
