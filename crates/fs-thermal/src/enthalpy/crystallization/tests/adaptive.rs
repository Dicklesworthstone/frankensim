//! G1 adiabatic recalescence and latent reheating; G4 transactional refusal.
use super::*;

fn first_order_law() -> CrystallizationLaw<'static> {
    CrystallizationLaw::new(
        &AMORPHOUS,
        &CRYSTALLINE,
        &RATE,
        CrystallizationConfig {
            identity: [0xAC; 32],
            avrami_exponent: 1.0,
            ..law().config()
        },
    )
    .unwrap()
}

fn adaptivity() -> CrystallizationAdaptivity {
    CrystallizationAdaptivity {
        maximum_dt_s: 20.0,
        minimum_dt_s: 1e-8,
        maximum_trials: 100_000,
        temperature_tolerance_k: 1e-5,
        crystalline_mass_fraction_tolerance: 1e-7,
        exposure_tolerance: 1e-6,
    }
}

#[test]
fn adaptive_history_resolves_analytic_recalescence_inside_one_thermal_step() {
    let law = first_order_law();
    let original =
        CrystallizationSimulation::new(thermal(&[380.0], 0.0, 0.0), &[material(&law)], || Ok(()))
            .unwrap();
    let selected = original
        .clone()
        .with_adaptive_phase_integration(MaterialId(3), adaptivity())
        .unwrap();
    let (next, audit) = selected
        .advance_to(20.0, controls(20.0), || Ok(()))
        .unwrap();

    // At fixed h=80000, T=380+100X. With n=1 and K=.002(500-T),
    // y=exp(z) satisfies y'=.08y+.16, y(0)=1. This reference solves
    // the coupled caloric/kinetic ODE independently of the midpoint kernel.
    let y = 3.0 * (0.08_f64 * 20.0).exp() - 2.0;
    let expected_fraction = 0.8 * (1.0 - 1.0 / y);
    let expected_temperature = 380.0 + 100.0 * expected_fraction;
    close(
        next.temperature_k(CellId(10)).unwrap(),
        expected_temperature,
        0.002,
    );
    close(
        next.states()[&CellId(10)].crystalline_mass_fraction(),
        expected_fraction,
        2e-5,
    );
    assert_eq!(next.energy_j(), original.energy_j());
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    assert_eq!(audit.thermal.stored_energy_change_j, 0.0);
    assert_eq!(audit.thermal.substeps, 1);
    assert!(audit.adaptive_phase.rejected_trials > 0);
    assert!(audit.adaptive_phase.accepted_substeps > 2);
    assert!(
        audit.adaptive_phase.maximum_temperature_discrepancy_k
            <= adaptivity().temperature_tolerance_k
    );
    assert!(
        audit
            .adaptive_phase
            .maximum_crystalline_mass_fraction_discrepancy
            <= adaptivity().crystalline_mass_fraction_tolerance
    );
    assert!(audit.adaptive_phase.maximum_exposure_discrepancy <= adaptivity().exposure_tolerance);
    assert_eq!(next.phase_adaptivity(MaterialId(3)), Some(adaptivity()));

    // A large unchanged outer step no longer determines kinetic accuracy.
    let fixed = original
        .advance_to(20.0, controls(20.0), || Ok(()))
        .unwrap()
        .0;
    let fixed_error = (fixed.temperature_k(CellId(10)).unwrap() - expected_temperature).abs();
    let adaptive_error = (next.temperature_k(CellId(10)).unwrap() - expected_temperature).abs();
    assert!(adaptive_error < 0.01 * fixed_error);
    assert_eq!(original.phase_adaptivity(MaterialId(3)), None);
    let replay = selected
        .advance_to(20.0, controls(20.0), || Ok(()))
        .unwrap();
    assert_eq!(next.cells(), replay.0.cells());
    assert_eq!(next.states(), replay.0.states());
    assert_eq!(audit, replay.1);
}

#[test]
fn adaptive_reheating_uses_the_same_heat_once_and_preserves_a_plateau_hold() {
    let law = first_order_law();
    let remelting = law.equivalent_age_remelting().unwrap();
    let selected =
        CrystallizationSimulation::new(thermal(&[380.0], 0.0, 0.0), &[material(&law)], || Ok(()))
            .unwrap()
            .with_equilibrium_remelting(MaterialId(3), remelting)
            .unwrap()
            .with_adaptive_phase_integration(MaterialId(3), adaptivity())
            .unwrap();
    let (mut aged, _) = selected
        .advance_to(20.0, controls(20.0), || Ok(()))
        .unwrap();
    let before = aged.states()[&CellId(10)].crystalline_mass_fraction();
    assert!(before > 0.3);
    let mut heater = boundary(0.0);
    heater.source_w = 100.0;
    aged.set_boundary(CellId(10), heater, || Ok(())).unwrap();
    let (mut melted, audit) = aged.advance_to(29.0, controls(9.0), || Ok(())).unwrap();

    // 0.01 kg receives 900 J: h=170000, T=500, X=(200000-h)/100000=.3.
    close(melted.cells()[&CellId(10)].specific_h_j_kg, 170_000.0, 1e-9);
    close(melted.temperature_k(CellId(10)).unwrap(), 500.0, 1e-10);
    close(
        melted.states()[&CellId(10)].crystalline_mass_fraction(),
        0.3,
        2e-12,
    );
    close(audit.thermal.external_heat_j, 900.0, 1e-10);
    close(audit.thermal.stored_energy_change_j, 900.0, 1e-10);
    close(audit.crystallized_mass_kg, 0.01 * (0.3 - before), 2e-14);
    assert_eq!(melted.mass_kg(), aged.mass_kg());
    assert_eq!(audit.thermal.substeps, 1);
    // The same plateau endpoint can be exact for one trial and its halves.
    // This endpoint error estimate does not claim a resolved melting event time.
    assert!(audit.adaptive_phase.accepted_substeps >= 2);

    melted
        .set_boundary(CellId(10), boundary(0.0), || Ok(()))
        .unwrap();
    let (held, hold_audit) = melted.advance_to(30.0, controls(1.0), || Ok(())).unwrap();
    assert_eq!(held.states(), melted.states());
    assert_eq!(hold_audit.crystallized_mass_kg, 0.0);
    assert_eq!(hold_audit.thermal.external_heat_j, 0.0);
    assert_eq!(held.remelting_model(MaterialId(3)), Some(remelting));
    assert_eq!(held.phase_adaptivity(MaterialId(3)), Some(adaptivity()));
}

#[test]
fn adaptive_selection_budget_and_late_cancellation_preserve_accepted_state() {
    let law = first_order_law();
    let original =
        CrystallizationSimulation::new(thermal(&[380.0], 0.0, 0.0), &[material(&law)], || Ok(()))
            .unwrap();
    assert!(
        original
            .clone()
            .with_adaptive_phase_integration(MaterialId(999), adaptivity())
            .is_err()
    );
    assert!(
        original
            .clone()
            .with_adaptive_phase_integration(
                MaterialId(3),
                CrystallizationAdaptivity {
                    maximum_trials: 0,
                    ..adaptivity()
                }
            )
            .is_err()
    );
    let selected = original
        .clone()
        .with_adaptive_phase_integration(MaterialId(3), adaptivity())
        .unwrap();
    assert!(
        selected
            .clone()
            .with_adaptive_phase_integration(MaterialId(3), adaptivity())
            .is_err()
    );

    let limited = original
        .clone()
        .with_adaptive_phase_integration(
            MaterialId(3),
            CrystallizationAdaptivity {
                maximum_trials: 1,
                ..adaptivity()
            },
        )
        .unwrap();
    assert!(matches!(
        limited.advance_to(20.0, controls(20.0), || Ok(())),
        Err(CrystallizationTimelineError::Material(
            CrystallizationError::NonConverged
        ))
    ));
    assert_eq!(limited.states(), original.states());
    assert_eq!(limited.cells(), original.cells());
    assert_eq!(limited.time_s(), original.time_s());

    let mut calls = 0;
    let refused = selected.advance_to(20.0, controls(20.0), || {
        calls += 1;
        if calls == 40 {
            Err(DepositionError::Budget("caller interruption"))
        } else {
            Ok(())
        }
    });
    assert_eq!(calls, 40);
    assert!(matches!(
        refused,
        Err(CrystallizationTimelineError::Thermal(
            DepositionError::Budget("caller interruption")
        ))
    ));
    assert_eq!(selected.states(), original.states());
    assert_eq!(selected.cells(), original.cells());
    assert_eq!(selected.time_s(), original.time_s());
    assert_eq!(selected.pending_events(), original.pending_events());
}
