//! Independent references for direct rank-one spring/contact equilibrium.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::normal_patch::{
    FiniteGapContactError, FiniteGapGrid, FiniteGapPlasticConfig, FiniteGapPlasticDrive,
    FiniteGapPlasticDriveError, FiniteGapPlasticError, FiniteGapPlasticLinearSolver,
    FiniteGapPlasticLoadControls, FiniteGapPlasticModel,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

const LINEAR_SOLVERS: [FiniteGapPlasticLinearSolver; 2] = [
    FiniteGapPlasticLinearSolver::Dense,
    FiniteGapPlasticLinearSolver::MatrixFree {
        maximum_iterations: 128,
        maximum_corrections: 4,
    },
];

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 107,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn config() -> FiniteGapPlasticConfig {
    let mut gap = vec![0.01; 25];
    gap[12] = 0.0;
    FiniteGapPlasticConfig {
        grid: FiniteGapGrid {
            cells_x: 5,
            cells_y: 5,
            cell_width_m: 1e-3,
            cell_depth_m: 1e-3,
        },
        undeformed_gap_m: gap,
        reduced_modulus_pa: 1e9,
        hardness_pa: vec![1e6; 25],
        maximum_active_set_iterations: 256,
        complementarity_tolerance_m: 1e-12,
        boundary_clearance_cells: 1,
    }
}

fn controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -0.02,
        maximum_approach_m: 0.02,
        absolute_force_tolerance_n: 1e-13,
        relative_force_tolerance: 1e-12,
        approach_tolerance_m: 1e-18,
        maximum_equilibrium_evaluations: 1,
    }
}

fn drive(position_m: f64) -> FiniteGapPlasticDrive {
    FiniteGapPlasticDrive {
        stiffness_n_per_m: 2e5,
        position_m,
    }
}

fn self_compliance(config: &FiniteGapPlasticConfig) -> f64 {
    4.0 * config.grid.cell_width_m * (1.0 + 2.0_f64.sqrt()).ln()
        / (core::f64::consts::PI * config.reduced_modulus_pa)
}

fn adjacent_compliance(config: &FiniteGapPlasticConfig) -> f64 {
    let s = config.grid.cell_width_m;
    let b = 0.5 * s;
    let primitive = |x: f64| x * (b / x).asinh() + b * (x / b).asinh();
    2.0 * (primitive(1.5 * s) - primitive(0.5 * s))
        / (core::f64::consts::PI * config.reduced_modulus_pa)
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-25),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

#[test]
fn one_equilibrium_ignores_unsafe_search_endpoints_and_keeps_contact_energy_physical() {
    with_cx(|cx, _| {
        for linear_solver in LINEAR_SOLVERS {
            let config = config();
            let c = self_compliance(&config);
            let a = config.grid.cell_width_m * config.grid.cell_depth_m;
            let model =
                FiniteGapPlasticModel::new_with_linear_solver(config, linear_solver, cx).unwrap();
            let before = model.initial_state();
            let saved = before.clone();
            // This caller bound is deliberately not a physically admissible contact
            // trial. It must not obstruct the actual interior spring equilibrium.
            assert!(matches!(
                model.advance(&before, controls().maximum_approach_m, cx),
                Err(FiniteGapPlasticError::Contact(
                    FiniteGapContactError::ContactReachedDomainBoundary { .. }
                ))
            ));
            let loading = drive(3e-6);
            let next = model
                .advance_with_drive(&before, loading, controls(), cx)
                .unwrap();
            let p = loading.position_m / (c + a / loading.stiffness_n_per_m);
            close(next.contact.candidate.pressure_pa()[12], p, 4e-12);
            close(next.contact.candidate.approach_m(), c * p, 4e-12);
            close(
                next.contact.receipt.normal.elastic_opening_m[12],
                c * p,
                4e-12,
            );
            close(
                next.contact.receipt.normal.reversible_energy_j,
                0.5 * a * c * p * p,
                4e-12,
            );
            close(
                next.drive.spring_energy_j,
                0.5 * (a * p) * (a * p) / loading.stiffness_n_per_m,
                4e-12,
            );
            assert_eq!(next.contact.receipt.plastic_work_j, 0.0);
            assert_eq!(next.contact.candidate.plastic_depth_m(), &[0.0; 25]);
            assert_eq!(next.drive.equilibrium_evaluations, 1);
            assert_eq!(next.drive.approach_bracket_m, [next.drive.approach_m; 2]);
            assert_eq!(
                next.drive.balance_residual_bracket_n,
                [next.drive.force_residual_n; 2]
            );
            assert_eq!(before, saved);
        }
    });
}

#[test]
fn unequal_hardness_caps_obey_coupled_equilibrium_and_original_plastic_work() {
    with_cx(|cx, _| {
        for linear_solver in LINEAR_SOLVERS {
            let mut config = config();
            config.undeformed_gap_m[13] = 2e-7;
            config.hardness_pa[12] = 2e5;
            let c = self_compliance(&config);
            let cross = adjacent_compliance(&config);
            let a = config.grid.cell_width_m * config.grid.cell_depth_m;
            let h = config.hardness_pa[12];
            let loading = drive(3e-6);
            let spring_compliance = a / loading.stiffness_n_per_m;
            // Cell 12 saturates; the adjacent, harder cell remains elastic. Its
            // equality includes both the half-space coupling and shared spring.
            let second_pressure = (loading.position_m
                - config.undeformed_gap_m[13]
                - (cross + spring_compliance) * h)
                / (c + spring_compliance);
            assert!(second_pressure > 0.0 && second_pressure < config.hardness_pa[13]);
            let approach = loading.position_m - spring_compliance * (h + second_pressure);
            let first_opening = c * h + cross * second_pressure;
            let second_opening = cross * h + c * second_pressure;
            let dent = approach - first_opening;
            assert!(dent > 0.0);
            let model =
                FiniteGapPlasticModel::new_with_linear_solver(config, linear_solver, cx).unwrap();
            let before = model.initial_state();
            let next = model
                .advance_with_drive(&before, loading, controls(), cx)
                .unwrap();
            assert_eq!(next.contact.candidate.pressure_pa()[12], h);
            close(
                next.contact.candidate.pressure_pa()[13],
                second_pressure,
                5e-12,
            );
            close(next.contact.candidate.approach_m(), approach, 5e-12);
            close(next.contact.candidate.plastic_depth_m()[12], dent, 5e-12);
            assert_eq!(next.contact.candidate.plastic_depth_m()[13], 0.0);
            close(
                next.contact.receipt.normal.elastic_opening_m[12],
                first_opening,
                5e-12,
            );
            close(
                next.contact.receipt.normal.elastic_opening_m[13],
                second_opening,
                5e-12,
            );
            close(next.contact.receipt.plastic_work_j, a * h * dent, 5e-12);
            close(
                next.contact.receipt.normal.reversible_energy_j,
                0.5 * a * (h * first_opening + second_pressure * second_opening),
                5e-12,
            );
            assert_eq!(next.contact.receipt.yielded_cells, 1);
            assert_eq!(next.drive.equilibrium_evaluations, 1);
            assert_eq!(next.contact.receipt.predecessor_identity, before.identity());
            // Independently prescribe the analytic approach through the original
            // contact owner, without any loading-spring compliance in that solve.
            let oracle = model.advance(&before, approach, cx).unwrap();
            for (actual, expected) in next
                .contact
                .candidate
                .pressure_pa()
                .iter()
                .zip(oracle.candidate.pressure_pa())
            {
                close(*actual, *expected, 6e-12);
            }
            close(
                next.contact.receipt.plastic_work_j,
                oracle.receipt.plastic_work_j,
                6e-12,
            );
            assert_eq!(before.plastic_depth_m(), &[0.0; 25]);
            assert_eq!(before.cumulative_plastic_work_j(), 0.0);
        }
    });
}

#[test]
fn bounds_admit_the_selected_point_without_changing_the_pressure_solution() {
    with_cx(|cx, _| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let before = model.initial_state();
        let loading = drive(3e-6);
        let wide = model
            .advance_with_drive(&before, loading, controls(), cx)
            .unwrap();
        let mut point = controls();
        point.minimum_approach_m = wide.drive.approach_m;
        point.maximum_approach_m = wide.drive.approach_m;
        point.approach_tolerance_m = 1.0;
        let exact = model
            .advance_with_drive(&before, loading, point, cx)
            .unwrap();
        assert_eq!(exact.contact, wide.contact);
        assert_eq!(exact.drive.equilibrium_evaluations, 1);
        assert_eq!(exact.drive.approach_bracket_m, [wide.drive.approach_m; 2]);
        point.minimum_approach_m = wide.drive.approach_m.next_up();
        point.maximum_approach_m = point.minimum_approach_m;
        match model.advance_with_drive(&before, loading, point, cx) {
            Err(FiniteGapPlasticDriveError::TargetOutsideApproachBounds {
                approach_bracket_m,
                balance_residual_bracket_n,
            }) => {
                assert_eq!(approach_bracket_m, [wide.drive.approach_m; 2]);
                assert_eq!(balance_residual_bracket_n, [wide.drive.force_residual_n; 2]);
            }
            other => panic!("expected actual selected point outside bounds, got {other:?}"),
        }
        assert_eq!(before, model.initial_state());
    });
}

#[test]
fn exact_detachment_keeps_the_unloaded_position_with_one_equilibrium() {
    with_cx(|cx, _| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let loading = drive(-1e-6);
        let mut point = controls();
        point.minimum_approach_m = loading.position_m;
        point.maximum_approach_m = loading.position_m;
        let next = model
            .advance_with_drive(&model.initial_state(), loading, point, cx)
            .unwrap();
        assert_eq!(next.contact.candidate.approach_m(), loading.position_m);
        assert_eq!(next.contact.candidate.pressure_pa(), &[0.0; 25]);
        assert_eq!(next.drive.drive_force_n, 0.0);
        assert_eq!(next.drive.spring_energy_j, 0.0);
        assert_eq!(next.drive.equilibrium_evaluations, 1);
    });
}

#[test]
fn direct_algebra_does_not_hide_an_unresolved_representable_force_residual() {
    with_cx(|cx, _| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let before = model.initial_state();
        let loading = FiniteGapPlasticDrive {
            stiffness_n_per_m: 3e5,
            position_m: core::f64::consts::PI * 1e-5,
        };
        let next = model
            .advance_with_drive(&before, loading, controls(), cx)
            .unwrap();
        // A saturated cell integrates to exactly 1 N. Reconstructing the spring
        // from the represented approach has a small, genuine nonzero residual.
        assert_eq!(next.drive.achieved_force_n, 1.0);
        let actual_spring =
            loading.stiffness_n_per_m * (loading.position_m - next.contact.candidate.approach_m());
        assert_eq!(next.drive.drive_force_n, actual_spring);
        assert_eq!(next.drive.force_residual_n, 1.0 - actual_spring);
        assert_ne!(next.drive.force_residual_n, 0.0);
        let mut strict = controls();
        strict.absolute_force_tolerance_n = 0.5 * next.drive.force_residual_n.abs();
        strict.relative_force_tolerance = 0.0;
        match model.advance_with_drive(&before, loading, strict, cx) {
            Err(FiniteGapPlasticDriveError::DidNotConverge {
                equilibrium_evaluations,
                approach_bracket_m,
                best_force_residual_n,
            }) => {
                assert_eq!(equilibrium_evaluations, 1);
                assert_eq!(approach_bracket_m, [next.drive.approach_m; 2]);
                assert_eq!(
                    best_force_residual_n,
                    Some(next.drive.force_residual_n.abs())
                );
            }
            other => panic!("expected strict physical force-balance refusal, got {other:?}"),
        }
        assert_eq!(before, model.initial_state());
    });
}

#[test]
fn one_equilibrium_still_honors_native_active_set_budget_and_cancellation() {
    with_cx(|cx, gate| {
        let mut bounded = config();
        bounded.maximum_active_set_iterations = 1;
        let model = FiniteGapPlasticModel::new(bounded, cx).unwrap();
        let before = model.initial_state();
        let saved = before.clone();
        assert!(matches!(
            model.advance_with_drive(&before, drive(1e-5), controls(), cx),
            Err(FiniteGapPlasticDriveError::Contact(
                FiniteGapPlasticError::Contact(FiniteGapContactError::ActiveSetDidNotConverge {
                    iterations: 1,
                    ..
                })
            ))
        ));
        assert_eq!(before, saved);
        gate.request();
        assert!(matches!(
            model.advance_with_drive(&before, drive(3e-6), controls(), cx),
            Err(FiniteGapPlasticDriveError::Contact(
                FiniteGapPlasticError::Contact(FiniteGapContactError::Cancelled)
            ))
        ));
        assert_eq!(before, saved);
    });
}
