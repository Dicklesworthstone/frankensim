//! Optional pressure algorithms retain original hardness, history, and bounds.

use super::*;
use fs_contact::normal_patch::{
    FiniteGapContactError, FiniteGapPlasticDrive, FiniteGapPlasticError,
    FiniteGapPlasticLinearSolver, FiniteGapPlasticLoadControls,
};

const MATRIX_FREE: FiniteGapPlasticLinearSolver = FiniteGapPlasticLinearSolver::MatrixFree {
    maximum_iterations: 128,
    maximum_corrections: 2,
};

fn rebuild(
    original: &SourcePlasticContact,
    envelope: PlasticContactEnvelope,
    selection: FiniteGapPlasticLinearSolver,
    cx: &Cx<'_>,
) -> Result<SourcePlasticContact, SourcePlasticContactError> {
    let config = original.model().config();
    SourcePlasticContact::new_with_linear_solver(
        original.normal_source().clone(),
        original.hardness_sources().clone(),
        config.grid,
        config.undeformed_gap_m.clone(),
        PlasticContactControls {
            maximum_active_set_iterations: config.maximum_active_set_iterations,
            complementarity_tolerance_m: config.complementarity_tolerance_m,
            boundary_clearance_cells: config.boundary_clearance_cells,
        },
        envelope,
        selection,
        cx,
    )
}

fn boundary_controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1.0e-4,
        maximum_approach_m: 1.0e-4,
        absolute_force_tolerance_n: 1.0e-14,
        relative_force_tolerance: 1.0e-12,
        approach_tolerance_m: 1.0e-18,
        maximum_equilibrium_evaluations: 128,
    }
}

#[test]
fn matrix_free_retains_original_sources_and_elastic_plastic_histories() {
    with_cx(|cx, _| {
        let (normal, hardness) = sources([1.0e6, 2.0e6], "linear-history", 0.0);
        let dense = build(normal, hardness, envelope(), cx).unwrap();
        let iterative = rebuild(&dense, envelope(), MATRIX_FREE, cx).unwrap();
        assert_eq!(
            dense.model().linear_solver(),
            FiniteGapPlasticLinearSolver::Dense
        );
        assert_eq!(iterative.model().linear_solver(), MATRIX_FREE);
        assert_eq!(iterative.normal_source(), dense.normal_source());
        assert_eq!(iterative.hardness_sources(), dense.hardness_sources());
        assert_eq!(
            rebuild(&dense, envelope(), FiniteGapPlasticLinearSolver::Dense, cx).unwrap(),
            dense
        );
        assert_ne!(iterative.model().identity(), dense.model().identity());
        assert_ne!(iterative.identity(), dense.identity());

        let mut expected = dense.initial_state();
        let mut actual = iterative.initial_state();
        for approach_m in [1.0e-6, 5.0e-5, -1.0e-5, 5.0e-5] {
            let reference = dense.advance(&expected, approach_m, cx).unwrap();
            let candidate = iterative.advance(&actual, approach_m, cx).unwrap();
            for (a, b) in [
                (
                    candidate.receipt.normal.normal_force_n,
                    reference.receipt.normal.normal_force_n,
                ),
                (
                    candidate.receipt.normal.reversible_energy_j,
                    reference.receipt.normal.reversible_energy_j,
                ),
                (
                    candidate.receipt.plastic_work_j,
                    reference.receipt.plastic_work_j,
                ),
                (
                    candidate.receipt.cumulative_plastic_work_j,
                    reference.receipt.cumulative_plastic_work_j,
                ),
            ] {
                // A repeated plastic endpoint may differ by roundoff-sized
                // work even when one backend returns exactly zero.
                assert!((a - b).abs() <= 2.0e-11 * a.abs().max(b.abs()) + 1.0e-18);
            }
            let got = candidate.candidate.contact_state();
            let want = reference.candidate.contact_state();
            for (a, b) in got.pressure_pa().iter().zip(want.pressure_pa()) {
                common::close(*a, *b);
            }
            for (a, b) in got.plastic_depth_m().iter().zip(want.plastic_depth_m()) {
                common::close(*a, *b);
            }
            assert_eq!(got.recession_m(), want.recession_m());
            assert_eq!(candidate.source_identity, iterative.identity());
            expected = reference.candidate;
            actual = candidate.candidate;
        }
        assert!(actual.contact_state().plastic_depth_m()[12] > 0.0);
    });
}

#[test]
fn numerical_selection_and_original_test_both_own_source_history() {
    with_cx(|cx, _| {
        let (normal, hardness) = sources([1.0e6, 2.0e6], "linear-owner", 0.0);
        let dense = build(normal, hardness, envelope(), cx).unwrap();
        let iterative = rebuild(&dense, envelope(), MATRIX_FREE, cx).unwrap();
        let limited = rebuild(
            &dense,
            envelope(),
            FiniteGapPlasticLinearSolver::MatrixFree {
                maximum_iterations: 64,
                maximum_corrections: 2,
            },
            cx,
        )
        .unwrap();
        let (normal, hardness) = sources([1.0e6, 2.0e6], "another-indentation-test", 0.0);
        let foreign = build(normal, hardness, envelope(), cx).unwrap();
        let foreign = rebuild(&foreign, envelope(), MATRIX_FREE, cx).unwrap();
        assert_eq!(foreign.model().identity(), iterative.model().identity());
        assert_ne!(foreign.identity(), iterative.identity());
        let before = iterative.initial_state();
        let saved = before.clone();
        for other in [&dense, &limited, &foreign] {
            assert!(matches!(
                other.advance(&before, 1.0e-6, cx),
                Err(SourcePlasticContactError::Invalid(_))
            ));
        }
        assert_eq!(before, saved);
        assert!(iterative.advance(&before, 1.0e-6, cx).is_ok());
    });
}

#[test]
fn matrix_free_retains_physical_envelope_and_native_cancellation() {
    with_cx(|cx, gate| {
        let (normal, hardness) = sources([1.0e6, 2.0e6], "linear-envelope", 0.0);
        let dense = build(normal, hardness, envelope(), cx).unwrap();
        for (bounds, quantity) in [
            (
                PlasticContactEnvelope {
                    maximum_plastic_depth_m: 1.0e-6,
                    ..envelope()
                },
                "permanent relative indentation [m]",
            ),
            (
                PlasticContactEnvelope {
                    half_space_depth_m: [1.0e-5; 2],
                    ..envelope()
                },
                "patch radius / homogeneous depth",
            ),
        ] {
            let law = rebuild(&dense, bounds, MATRIX_FREE, cx).unwrap();
            let before = law.initial_state();
            let saved = before.clone();
            assert!(matches!(
                law.advance(&before, 5.0e-5, cx),
                Err(SourcePlasticContactError::OutsideApplicability {
                    quantity: rejected,
                    ..
                }) if rejected == quantity
            ));
            assert_eq!(before, saved);
        }
        let law = rebuild(&dense, envelope(), MATRIX_FREE, cx).unwrap();
        let before = law.initial_state();
        let saved = before.clone();
        gate.request();
        assert!(matches!(
            law.advance(&before, 1.0e-6, cx),
            Err(SourcePlasticContactError::Contact(
                FiniteGapPlasticError::Contact(FiniteGapContactError::Cancelled)
            ))
        ));
        assert_eq!(before, saved);
    });
}

#[test]
fn every_thermal_boundary_preserves_the_pressure_solver_and_original_claims() {
    with_cx(|cx, _| {
        let tag = "linear-temperature";
        let cards = [material("A", 1.0e6, tag).0, material("B", 2.0e6, tag).0];
        let (normal, hardness) = sources([1.0e6, 2.0e6], tag, 0.0);
        let dense = build(normal, hardness, envelope(), cx).unwrap();
        for selection in [
            MATRIX_FREE,
            FiniteGapPlasticLinearSolver::Fft {
                maximum_iterations: 128,
                maximum_corrections: 2,
            },
        ] {
            let law = rebuild(&dense, envelope(), selection, cx).unwrap();
            let dented = law.advance(&law.initial_state(), 5.0e-5, cx).unwrap();
            let force = 0.005;
            let elastic = law
                .advance_at_load(&dented.candidate, force, boundary_controls(), cx)
                .unwrap();
            let before = elastic.candidate;
            let saved = before.clone();
            let drive = FiniteGapPlasticDrive {
                stiffness_n_per_m: 1000.0,
                position_m: before.contact_state().approach_m() + force / 1000.0,
            };
            // Constant-valued claims isolate reconstruction from changes in physics;
            // the independently changed query coordinates still create new owners.
            let hot = [313.15, 303.15];
            let held = law
                .refresh_at_temperatures(&before, [&cards[0], &cards[1]], hot, 308.15, cx)
                .unwrap();
            let loaded = law
                .refresh_at_temperatures_at_load(
                    &before,
                    [&cards[0], &cards[1]],
                    hot,
                    308.15,
                    force,
                    boundary_controls(),
                    cx,
                )
                .unwrap();
            let driven = law
                .refresh_at_temperatures_with_drive(
                    &before,
                    [&cards[0], &cards[1]],
                    hot,
                    308.15,
                    drive,
                    boundary_controls(),
                    cx,
                )
                .unwrap();
            for (source, candidate) in [
                (held.source(), held.candidate()),
                (loaded.source(), loaded.candidate()),
                (driven.source(), driven.candidate()),
            ] {
                assert_eq!(source.model().linear_solver(), selection);
                assert_ne!(source.identity(), law.identity());
                assert_eq!(
                    candidate.contact_state().plastic_depth_m(),
                    before.contact_state().plastic_depth_m()
                );
                assert_eq!(
                    candidate.contact_state().recession_m(),
                    before.contact_state().recession_m()
                );
                for side in 0..2 {
                    let original = &law.hardness_sources()[side];
                    let refreshed = &source.hardness_sources()[side];
                    assert_eq!(refreshed.key(), original.key());
                    assert_eq!(
                        refreshed.resolved().card_identity(),
                        original.resolved().card_identity()
                    );
                    assert_eq!(refreshed.resolved().state_coordinate("T"), Some(hot[side]));
                    assert_eq!(
                        refreshed.property().answer().receipt.selected,
                        original.property().answer().receipt.selected
                    );
                }
                assert!(
                    source
                        .advance(candidate, candidate.contact_state().approach_m(), cx)
                        .is_ok()
                );
            }
            assert_eq!(before, saved);
        }
    });
}
