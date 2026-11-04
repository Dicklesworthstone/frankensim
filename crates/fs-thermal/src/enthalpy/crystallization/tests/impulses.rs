//! Zero-time calorimetry, original history and atomic batch refusal.
use super::*;

fn impulse_controls() -> CrystallizationImpulseControls {
    let original = controls(0.01);
    CrystallizationImpulseControls {
        caloric: original.caloric,
        cell_energy_tolerance_j: original.cell_energy_tolerance_j,
        total_energy_tolerance_j: original.total_energy_tolerance_j,
    }
}

#[test]
fn simultaneous_heat_and_cooling_preserve_time_and_resume_real_contact_history() {
    let law = law();
    let bath = finite_bath();
    let mut events: Vec<_> = bath
        .cells()
        .values()
        .cloned()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    events.extend(bath.links().values().copied().map(|link| TimedEvent {
        time_s: 0.0,
        event: Event::Contact(link),
    }));
    events.push(TimedEvent {
        time_s: 1.0,
        event: Event::Separation(1),
    });
    let start = CrystallizationSimulation::new(
        activate(
            events,
            &[(MaterialId(3), &AMORPHOUS), (MaterialId(4), &MOLD)],
        ),
        &[material(&law)],
        || Ok(()),
    )
    .unwrap();
    let heat = BTreeMap::from([(CellId(10), 100.0), (CellId(11), -50.0)]);
    let (next, audit) = start
        .apply_heat_impulses(&heat, impulse_controls(), || Ok(()))
        .unwrap();

    // Independent calorimetry: 10g at Cp=1000 gains 10K; 50g at Cp=500 loses 2K.
    close(next.temperature_k(CellId(10)).unwrap(), 390.0, 1e-9);
    close(next.temperature_k(CellId(11)).unwrap(), 298.0, 1e-12);
    assert_eq!(next.states()[&CellId(10)].exposure(), 0.0);
    assert_eq!(next.states()[&CellId(10)].crystalline_mass_fraction(), 0.0);
    assert!(!next.states().contains_key(&CellId(11)));
    assert!(std::ptr::eq(
        next.material_law(MaterialId(3)).unwrap(),
        &law
    ));
    assert!(next.material_law(MaterialId(4)).is_none());
    assert_eq!(next.time_s(), start.time_s());
    assert_eq!(next.mass_kg(), start.mass_kg());
    assert_eq!(next.links(), start.links());
    assert_eq!(next.pending_events(), start.pending_events());
    assert_eq!(next.pending_events().len(), 1);
    for (&id, cell) in next.cells() {
        let mut expected = start.cells()[&id].clone();
        expected.specific_h_j_kg = cell.specific_h_j_kg;
        assert_eq!(*cell, expected);
    }
    close(next.energy_j() - start.energy_j(), 50.0, 1e-10);
    assert_eq!(audit.thermal.external_heat_j, 50.0);
    close(audit.thermal.stored_energy_change_j, 50.0, 1e-10);
    assert_eq!(audit.thermal.incoming_mass_kg, 0.0);
    assert_eq!(audit.thermal.incoming_energy_j, 0.0);
    assert_eq!((audit.thermal.substeps, audit.thermal.events), (0, 0));
    assert_eq!(audit.crystallized_mass_kg, 0.0);
    assert_eq!(
        audit.adaptive_phase,
        CrystallizationAdaptiveAudit::default()
    );
    let replay = start
        .apply_heat_impulses(&heat, impulse_controls(), || Ok(()))
        .unwrap();
    assert_eq!(next.states(), replay.0.states());
    assert_eq!(next.cells(), replay.0.cells());
    assert_eq!(audit, replay.1);

    let (held, hold_audit) = next
        .apply_heat_impulses(
            &BTreeMap::from([(CellId(10), 0.0), (CellId(11), 0.0)]),
            impulse_controls(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(held.states(), next.states());
    assert_eq!(held.cells(), next.cells());
    assert_eq!(hold_audit.thermal.external_heat_j, 0.0);
    let (evolved, _) = held.advance_to(0.05, controls(0.01), || Ok(())).unwrap();
    assert!(evolved.states()[&CellId(10)].exposure() > 0.0);
    assert!(evolved.temperature_k(CellId(11)).unwrap() > 298.0);
    close(evolved.energy_j(), next.energy_j(), 1e-8);
    assert_eq!(evolved.pending_events(), start.pending_events());
    assert_eq!(start.states()[&CellId(10)].exposure(), 0.0);
    assert_eq!(start.cells()[&CellId(10)].specific_h_j_kg, 80_000.0);
}

#[test]
fn an_instantaneous_heat_credit_remelts_once_without_growing_during_recooling() {
    let law = law();
    let start =
        CrystallizationSimulation::new(thermal(&[380.0], 0.0, 0.0), &[material(&law)], || Ok(()))
            .unwrap();
    let (aged, _) = start.advance_to(20.0, controls(0.02), || Ok(())).unwrap();
    let old = aged.states()[&CellId(10)];
    assert!(old.crystalline_mass_fraction() > 0.3);
    let heat = BTreeMap::from([(CellId(10), 900.0)]);
    assert!(
        aged.apply_heat_impulses(&heat, impulse_controls(), || Ok(()))
            .is_err()
    );
    let selected = aged
        .clone()
        .with_equilibrium_remelting(MaterialId(3), law.equivalent_age_remelting().unwrap())
        .unwrap();
    let (melted, audit) = selected
        .apply_heat_impulses(&heat, impulse_controls(), || Ok(()))
        .unwrap();
    let state = melted.states()[&CellId(10)];
    // h=80000+900/.01=170000; at Tm, X=(200000-h)/100000=.3.
    close(state.specific_enthalpy_j_kg(), 170_000.0, 1e-9);
    close(state.temperature_k(), 500.0, 1e-10);
    close(state.crystalline_mass_fraction(), 0.3, 2e-12);
    assert!(state.exposure() < old.exposure());
    assert_eq!(melted.time_s(), 20.0);
    assert_eq!(melted.mass_kg(), aged.mass_kg());
    assert_eq!(audit.thermal.external_heat_j, 900.0);
    close(audit.thermal.stored_energy_change_j, 900.0, 1e-10);
    close(
        audit.crystallized_mass_kg,
        0.01 * (0.3 - old.crystalline_mass_fraction()),
        2e-14,
    );
    assert_eq!((audit.thermal.substeps, audit.thermal.events), (0, 0));
    let (cooled, cooling) = melted
        .apply_heat_impulses(
            &BTreeMap::from([(CellId(10), -900.0)]),
            impulse_controls(),
            || Ok(()),
        )
        .unwrap();
    let cool = cooled.states()[&CellId(10)];
    assert_eq!(cool.exposure(), state.exposure());
    assert_eq!(
        cool.crystalline_mass_fraction(),
        state.crystalline_mass_fraction()
    );
    close(cool.temperature_k(), 410.0, 1e-9);
    assert_eq!(cooling.crystallized_mass_kg, 0.0);
    assert_eq!(cooled.time_s(), 20.0);
    close(cooled.energy_j(), aged.energy_j(), 1e-10);
    let (regrown, _) = cooled.advance_to(20.1, controls(0.01), || Ok(())).unwrap();
    assert!(
        regrown.states()[&CellId(10)].crystalline_mass_fraction()
            > cool.crystalline_mass_fraction()
    );
    assert_eq!(aged.states()[&CellId(10)], old);
}

#[test]
fn late_refusal_cancellation_roundoff_and_due_events_publish_no_partial_batch() {
    let law = law();
    let start =
        CrystallizationSimulation::new(finite_bath(), &[material(&law)], || Ok(())).unwrap();
    let saved_cells = start.cells().clone();
    let saved_states = start.states().clone();
    let valid = BTreeMap::from([(CellId(10), 100.0), (CellId(11), -50.0)]);
    let mut successful_checkpoints = 0;
    start
        .apply_heat_impulses(&valid, impulse_controls(), || {
            successful_checkpoints += 1;
            Ok(())
        })
        .unwrap();
    let mut calls = 0;
    let cancelled = start.apply_heat_impulses(&valid, impulse_controls(), || {
        calls += 1;
        if calls == successful_checkpoints {
            Err(DepositionError::Budget("caller interruption"))
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        cancelled,
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::Budget("caller interruption")
        ))
    ));
    for invalid in [
        BTreeMap::from([(CellId(10), 100.0), (CellId(11), 1e9)]),
        BTreeMap::from([(CellId(10), 100.0), (CellId(99), 0.0)]),
        BTreeMap::from([(CellId(10), 100.0), (CellId(11), f64::NAN)]),
    ] {
        assert!(
            start
                .apply_heat_impulses(&invalid, impulse_controls(), || Ok(()))
                .is_err()
        );
    }
    let mut invalid_controls = impulse_controls();
    invalid_controls.caloric.maximum_iterations = 0;
    assert!(
        start
            .apply_heat_impulses(&BTreeMap::new(), invalid_controls, || Ok(()))
            .is_err()
    );

    // q/m is below one representable h increment; a rounded zero delta must
    // not satisfy a tighter energy budget merely because the state stayed put.
    let mut strict = impulse_controls();
    strict.cell_energy_tolerance_j = 1e-22;
    strict.total_energy_tolerance_j = 1e-22;
    assert_eq!(80_000.0_f64 + 1e-20 / 0.01, 80_000.0);
    assert!(matches!(
        start.apply_heat_impulses(&BTreeMap::from([(CellId(10), 1e-20)]), strict, || Ok(()),),
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::EnergyBalance(_)
        ))
    ));
    assert_eq!(start.cells(), &saved_cells);
    assert_eq!(start.states(), &saved_states);
    assert_eq!(start.time_s(), 0.0);

    let due = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(saved_cells[&CellId(10)].clone()),
        }],
        &[(MaterialId(3), &AMORPHOUS as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 1,
        },
    )
    .unwrap();
    let due = CrystallizationSimulation::new(due, &[material(&law)], || Ok(())).unwrap();
    assert!(matches!(
        due.apply_heat_impulses(&BTreeMap::new(), impulse_controls(), || Ok(())),
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::InvalidInput(_)
        ))
    ));
    assert!(due.cells().is_empty());
    assert_eq!(due.pending_events().len(), 1);
}
