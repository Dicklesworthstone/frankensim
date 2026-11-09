use super::*;
use crate::enthalpy::pyrolysis::implicit::adaptive::{AdaptiveControls, AdaptiveError};

fn accuracy(mass: f64) -> AdaptiveControls {
    AdaptiveControls {
        initial_dt_s: 1.0,
        minimum_dt_s: 1e-8,
        maximum_comparisons: 2000,
        mass_absolute_tolerance_kg: mass,
        relative_tolerance: 0.0,
        ..Default::default()
    }
}

#[test]
fn isothermal_reaction_accuracy_resolves_separate_virgin_char_and_actual_outflow() {
    let gas = |t| VIRGIN.specific_enthalpy_j_kg(t);
    for yield_char in [0.0, 0.4] {
        let mut old = setup(
            &VIRGIN,
            &gas,
            &[650.0],
            &[0.01],
            0.0,
            0.0,
            Some(&VIRGIN),
            yield_char,
            3.0,
        );
        old.thermal.budget.max_substeps = 4000;
        let exact_virgin = 0.01 * (-3.0_f64).exp();
        let (fixed, fixed_report) = advance(&old, 1.0, 1.0);
        assert!(fixed_report.pyrolysis.absolute_energy_residual_j < 1e-8);
        let fixed_error = (fixed.inventories()[&CellId(1)].virgin_mass_kg - exact_virgin).abs();
        let mut errors = [0.0; 2];
        let mut workspace = old.implicit_workspace().unwrap();
        for (rung, mass_tolerance) in [1e-5, 1e-7].into_iter().enumerate() {
            let (next, report) = old
                .advance_adaptive_implicit_to(
                    1.0,
                    PyrolysisControls {
                        maximum_dt_s: 1.0,
                        ..Default::default()
                    },
                    ImplicitSolve::default(),
                    accuracy(mass_tolerance),
                    &mut workspace,
                    || Ok(()),
                )
                .unwrap();
            let inventory = next.inventories()[&CellId(1)];
            errors[rung] = (inventory.virgin_mass_kg - exact_virgin).abs();
            close(
                inventory.char_mass_kg,
                yield_char * (0.01 - inventory.virgin_mass_kg),
                2e-16,
            );
            close(next.temperature_k(CellId(1)).unwrap(), 650.0, 1e-9);
            close(
                next.mass_kg() + report.physical.pyrolysis.volatile_mass_kg,
                old.mass_kg(),
                1e-15,
            );
            close(
                next.energy_j() + report.physical.pyrolysis.volatile_enthalpy_j,
                old.energy_j(),
                1e-8,
            );
            assert!(report.accuracy_rejections > 0);
            assert!(report.maximum_accepted_normalized_error <= 1.0);
            assert_eq!(
                report.physical.pyrolysis.thermal.substeps,
                2 * report.accepted_comparisons
            );
            assert_eq!(report.work.implicit_intervals, 3 * report.comparisons);
            assert_eq!(
                report.comparisons,
                report.accepted_comparisons + report.accuracy_rejections
            );
            assert!(report.work.local_cell_solves > report.physical.local_cell_solves);
        }
        assert!(
            errors[1] < 0.2 * errors[0],
            "yield={yield_char}: {errors:?}"
        );
        assert!(
            errors[1] < 0.01 * fixed_error && errors[1] < 1e-5,
            "{errors:?}"
        );
    }
}

#[test]
fn discarded_trials_emit_nothing_and_fine_receipt_matches_native_steps() {
    let gas = |t| VIRGIN.specific_enthalpy_j_kg(t);
    let old = setup(&VIRGIN, &gas, &[650.0], &[0.01], 0.0, 0.0, None, 0.0, 1.0);
    let (expected, native) = advance(&old, 1.0, 0.25);
    let mut workspace = old.implicit_workspace().unwrap();
    let (next, report) = old
        .advance_adaptive_implicit_to(
            1.0,
            PyrolysisControls {
                maximum_dt_s: 1.0,
                ..Default::default()
            },
            ImplicitSolve::default(),
            accuracy(3e-4),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(report.accuracy_rejections, 1);
    assert_eq!(report.accepted_comparisons, 2);
    assert_eq!(report.physical.pyrolysis.thermal.substeps, 4);
    assert_eq!(report.work.implicit_intervals, 9);
    assert_eq!(next.cells(), expected.cells());
    assert_eq!(next.inventories(), expected.inventories());
    close(next.mass_kg(), 0.01 / 1.25_f64.powi(4), 1e-16);
    close(
        report.physical.pyrolysis.volatile_mass_kg,
        native.pyrolysis.volatile_mass_kg,
        1e-16,
    );
    close(
        report.physical.pyrolysis.volatile_enthalpy_j,
        native.pyrolysis.volatile_enthalpy_j,
        1e-10,
    );
    assert_eq!(
        report.physical.nonlinear_iterations,
        native.nonlinear_iterations
    );
    assert_eq!(report.physical.pyrolysis.thermal.external_heat_j, 0.0);
    assert_eq!(report.physical.pyrolysis.releases.len(), 1);
}

#[test]
fn common_formation_reference_does_not_inflate_adaptive_energy_allowance() {
    let run = |offset| {
        let virgin = Caloric { cp: 2000.0, offset };
        let char_law = Caloric {
            cp: 1200.0,
            offset: offset + 400_000.0,
        };
        let gas = |t| Ok(offset + 2000.0 * (t - 300.0) + 200_000.0);
        let old = setup(
            &virgin,
            &gas,
            &[700.0, 500.0],
            &[0.01, 0.02],
            20.0,
            2.0,
            Some(&char_law),
            0.2,
            1.7,
        );
        let mut workspace = old.implicit_workspace().unwrap();
        let (next, report) = old
            .advance_adaptive_implicit_to(
                0.2,
                PyrolysisControls {
                    maximum_dt_s: 0.2,
                    cell_energy_tolerance_j: 1e-7,
                    total_energy_tolerance_j: 1e-5,
                    ..Default::default()
                },
                ImplicitSolve::default(),
                AdaptiveControls {
                    initial_dt_s: 0.2,
                    storage_absolute_tolerance_j_kg: 200.0,
                    temperature_absolute_tolerance_k: 0.1,
                    mass_absolute_tolerance_kg: 1e-6,
                    volatile_energy_absolute_tolerance_j: 0.1,
                    relative_tolerance: 0.0,
                    ..Default::default()
                },
                &mut workspace,
                || Ok(()),
            )
            .unwrap();
        (
            [
                next.temperature_k(CellId(1)).unwrap(),
                next.temperature_k(CellId(2)).unwrap(),
            ],
            next.mass_kg(),
            report,
        )
    };
    let original = run(0.0);
    let shifted = run(1e8);
    for i in 0..2 {
        close(original.0[i], shifted.0[i], 1e-6);
    }
    close(original.1, shifted.1, 1e-13);
    assert_eq!(original.2.comparisons, shifted.2.comparisons);
    assert_eq!(
        original.2.accepted_comparisons,
        shifted.2.accepted_comparisons
    );
    close(
        original.2.physical.pyrolysis.volatile_mass_kg,
        shifted.2.physical.pyrolysis.volatile_mass_kg,
        1e-13,
    );
    close(
        shifted.2.physical.pyrolysis.volatile_enthalpy_j
            - 1e8 * shifted.2.physical.pyrolysis.volatile_mass_kg,
        original.2.physical.pyrolysis.volatile_enthalpy_j,
        1e-6,
    );
}

#[test]
fn accuracy_comparison_and_original_substep_limits_refuse_the_whole_advance() {
    let gas = |t| VIRGIN.specific_enthalpy_j_kg(t);
    let mut old = setup(&VIRGIN, &gas, &[650.0], &[0.01], 0.0, 0.0, None, 0.0, 1.0);
    let cells = old.cells().clone();
    let inventories = old.inventories().clone();
    let controls = PyrolysisControls {
        maximum_dt_s: 1.0,
        ..Default::default()
    };
    let mut workspace = old.implicit_workspace().unwrap();
    for limit in [
        AdaptiveControls {
            minimum_dt_s: 1.0,
            ..accuracy(1e-12)
        },
        AdaptiveControls {
            maximum_reductions: 0,
            ..accuracy(1e-12)
        },
    ] {
        assert!(matches!(
            old.advance_adaptive_implicit_to(
                1.0,
                controls,
                ImplicitSolve::default(),
                limit,
                &mut workspace,
                || Ok(())
            ),
            Err(AdaptiveError::TimeAccuracy { .. })
        ));
    }
    assert_eq!(
        old.advance_adaptive_implicit_to(
            1.0,
            controls,
            ImplicitSolve::default(),
            AdaptiveControls {
                maximum_comparisons: 1,
                ..accuracy(3e-4)
            },
            &mut workspace,
            || Ok(())
        )
        .err(),
        Some(AdaptiveError::ComparisonBudget),
    );
    old.thermal.budget.max_substeps = 3;
    assert!(matches!(
        old.advance_adaptive_implicit_to(
            1.0,
            controls,
            ImplicitSolve::default(),
            accuracy(3e-4),
            &mut workspace,
            || Ok(())
        ),
        Err(AdaptiveError::Implicit(ImplicitError::Thermal(
            DepositionError::Budget(_)
        )))
    ));
    assert_eq!(old.cells(), &cells);
    assert_eq!(old.inventories(), &inventories);
    assert_eq!(old.time_s(), 0.0);
}

#[test]
fn late_cancellation_and_native_domain_refusal_preserve_state_and_workspace_reuse() {
    let gas = |t| VIRGIN.specific_enthalpy_j_kg(t);
    let old = setup(&VIRGIN, &gas, &[650.0], &[0.01], 0.0, 0.0, None, 0.0, 1.0);
    let controls = PyrolysisControls {
        maximum_dt_s: 1.0,
        ..Default::default()
    };
    let mut workspace = old.implicit_workspace().unwrap();
    let mut complete_polls = 0;
    let (expected, expected_report) = old
        .advance_adaptive_implicit_to(
            1.0,
            controls,
            ImplicitSolve::default(),
            accuracy(3e-4),
            &mut workspace,
            || {
                complete_polls += 1;
                Ok(())
            },
        )
        .unwrap();
    let mut polls = 0;
    let error = old
        .advance_adaptive_implicit_to(
            1.0,
            controls,
            ImplicitSolve::default(),
            accuracy(3e-4),
            &mut workspace,
            || {
                polls += 1;
                if polls == complete_polls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        )
        .err();
    assert_eq!(
        error,
        Some(AdaptiveError::Implicit(ImplicitError::Cancelled))
    );
    assert_eq!(old.time_s(), 0.0);
    close(old.mass_kg(), 0.01, 0.0);
    let (replayed, report) = old
        .advance_adaptive_implicit_to(
            1.0,
            controls,
            ImplicitSolve::default(),
            accuracy(3e-4),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(replayed.cells(), expected.cells());
    assert_eq!(replayed.inventories(), expected.inventories());
    assert_eq!(report, expected_report);

    let heated = setup(&VIRGIN, &gas, &[650.0], &[0.01], 1e6, 0.0, None, 0.0, 0.0);
    let native_error = heated
        .advance_implicit_to(
            1.0,
            controls,
            ImplicitSolve::default(),
            &mut workspace,
            || Ok(()),
        )
        .err()
        .unwrap();
    assert_eq!(
        heated
            .advance_adaptive_implicit_to(
                1.0,
                controls,
                ImplicitSolve::default(),
                accuracy(3e-4),
                &mut workspace,
                || Ok(())
            )
            .err(),
        Some(AdaptiveError::Implicit(native_error)),
    );
    assert_eq!(heated.time_s(), 0.0);
    close(heated.temperature_k(CellId(1)).unwrap(), 650.0, 0.0);
}
