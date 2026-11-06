//! Manufactured and analytical checks of the optional bounded pressure solve.
//! These check the numerical model; they are not material calibration data.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::normal_patch::{
    FiniteGapContactError, FiniteGapGrid, FiniteGapPlasticConfig, FiniteGapPlasticError,
    FiniteGapPlasticLinearSolver, FiniteGapPlasticModel, MAX_FINITE_GAP_PLASTIC_FREE_CELLS,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

const ITERATIVE: FiniteGapPlasticLinearSolver = FiniteGapPlasticLinearSolver::MatrixFree {
    maximum_iterations: 256,
    maximum_corrections: 2,
};

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 37,
                kernel_id: 181,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

// Independent asinh antiderivative of the rectangular 1/r integral. The
// implementation uses a logarithmic primitive and its own signed-offset cache.
fn rectangle_compliance(x: f64, y: f64, side: f64, modulus: f64) -> f64 {
    fn primitive(x: f64, y: f64) -> f64 {
        let a = if x == 0.0 {
            0.0
        } else {
            x * (y / x.abs()).asinh()
        };
        let b = if y == 0.0 {
            0.0
        } else {
            y * (x / y.abs()).asinh()
        };
        a + b
    }
    let h = 0.5 * side;
    (primitive(x + h, y + h) - primitive(x - h, y + h) - primitive(x + h, y - h)
        + primitive(x - h, y - h))
        / (core::f64::consts::PI * modulus)
}

fn small_config() -> FiniteGapPlasticConfig {
    let mut gap = vec![0.01; 25];
    gap[12] = 0.0;
    gap[13] = 1.0e-6;
    let mut hardness = vec![1.3e6; 25];
    hardness[12] = 0.4e6;
    FiniteGapPlasticConfig {
        grid: FiniteGapGrid {
            cells_x: 5,
            cells_y: 5,
            cell_width_m: 1.0e-3,
            cell_depth_m: 1.0e-3,
        },
        undeformed_gap_m: gap,
        reduced_modulus_pa: 1.0e9,
        hardness_pa: hardness,
        maximum_active_set_iterations: 64,
        complementarity_tolerance_m: 1.0e-12,
        boundary_clearance_cells: 1,
    }
}

#[test]
fn manufactured_patch_exceeds_dense_limit_with_known_pressure_and_force() {
    with_cx(|cx, _| {
        let n = 49_usize;
        let side = 1.0e-3;
        let modulus = 1.0e9;
        let pressure = 1.0e5;
        let active: Vec<usize> = (0..n * n)
            .filter(|i| i % n > 0 && i % n < n - 1 && i / n > 0 && i / n < n - 1)
            .collect();
        assert_eq!(active.len(), 2209);
        assert!(active.len() > MAX_FINITE_GAP_PLASTIC_FREE_CELLS);
        let kernel: Vec<f64> = (0..n * n)
            .map(|i| {
                rectangle_compliance((i % n) as f64 * side, (i / n) as f64 * side, side, modulus)
            })
            .collect();
        let mut gap = vec![1.0; n * n];
        let mut approach = 0.0_f64;
        for &i in &active {
            let opening: f64 = active
                .iter()
                .map(|&j| {
                    let dx = (i % n).abs_diff(j % n);
                    let dy = (i / n).abs_diff(j / n);
                    kernel[dy * n + dx] * pressure
                })
                .sum();
            gap[i] = -opening;
            approach = approach.max(opening);
        }
        let config = FiniteGapPlasticConfig {
            grid: FiniteGapGrid {
                cells_x: n,
                cells_y: n,
                cell_width_m: side,
                cell_depth_m: side,
            },
            undeformed_gap_m: gap,
            reduced_modulus_pa: modulus,
            hardness_pa: vec![1.0e9; n * n],
            maximum_active_set_iterations: 32,
            complementarity_tolerance_m: 1.0e-11,
            boundary_clearance_cells: 1,
        };
        let dense = FiniteGapPlasticModel::new(config.clone(), cx).unwrap();
        assert!(matches!(
            dense.advance(&dense.initial_state(), approach, cx),
            Err(FiniteGapPlasticError::FreeSetTooLarge {
                free_cells: 2209,
                maximum: 2048
            })
        ));
        let model =
            FiniteGapPlasticModel::new_with_linear_solver(config.clone(), ITERATIVE, cx).unwrap();
        let before = model.initial_state();
        let saved = before.clone();
        let step = model.advance(&before, approach, cx).unwrap();
        assert_eq!(before, saved);
        assert_eq!(step.receipt.normal.active_cells, active.len());
        assert_eq!(step.receipt.yielded_cells, 0);
        assert_eq!(step.receipt.plastic_work_j, 0.0);
        assert!(step.candidate.plastic_depth_m().iter().all(|z| *z == 0.0));
        assert!(step.receipt.kkt_residual_m <= config.complementarity_tolerance_m);
        assert!(
            step.receipt.normal.complementarity_residual_m <= config.complementarity_tolerance_m
        );
        close(
            step.receipt.normal.normal_force_n,
            active.len() as f64 * pressure * side * side,
            5.0e-6,
            0.0,
        );
        for (i, &actual) in step.candidate.pressure_pa().iter().enumerate() {
            if i % n == 0 || i % n == n - 1 || i / n == 0 || i / n == n - 1 {
                assert_eq!(actual, 0.0);
            } else {
                close(actual, pressure, 5.0e-5, 0.0);
                let independent_opening = -config.undeformed_gap_m[i];
                close(
                    step.receipt.normal.elastic_opening_m[i],
                    independent_opening,
                    0.0,
                    config.complementarity_tolerance_m,
                );
            }
        }
    });
}

#[test]
fn coupled_elastic_and_capped_cells_match_independent_equilibrium_and_work() {
    with_cx(|cx, _| {
        let config = small_config();
        let self_compliance = rectangle_compliance(0.0, 0.0, 1.0e-3, 1.0e9);
        let cross_compliance = rectangle_compliance(1.0e-3, 0.0, 1.0e-3, 1.0e9);
        let p_a = config.hardness_pa[12];
        let p_b = (1.0e-6 - cross_compliance * p_a) / self_compliance;
        let dent = 2.0e-6 - self_compliance * p_a - cross_compliance * p_b;
        assert!(p_b > 0.0 && p_b < config.hardness_pa[13]);
        assert!(dent > 0.0);
        for selection in [FiniteGapPlasticLinearSolver::Dense, ITERATIVE] {
            let model =
                FiniteGapPlasticModel::new_with_linear_solver(config.clone(), selection, cx)
                    .unwrap();
            let before = model.initial_state();
            let step = model.advance(&before, 2.0e-6, cx).unwrap();
            assert_eq!(step.candidate.pressure_pa()[12], p_a);
            close(step.candidate.pressure_pa()[13], p_b, 1.0e-10, 0.0);
            close(step.candidate.plastic_depth_m()[12], dent, 1.0e-10, 0.0);
            assert_eq!(step.candidate.plastic_depth_m()[13], 0.0);
            close(
                step.receipt.plastic_work_j,
                p_a * dent * 1.0e-6,
                1.0e-10,
                0.0,
            );
            assert_eq!(step, model.advance(&before, 2.0e-6, cx).unwrap());
            let unloaded = model.advance(&step.candidate, -1.0e-6, cx).unwrap();
            assert_eq!(unloaded.receipt.normal.normal_force_n, 0.0);
            assert_eq!(
                unloaded.candidate.plastic_depth_m(),
                step.candidate.plastic_depth_m()
            );
            let mut recession = vec![0.0; 25];
            recession[12] = 2.0e-7;
            let worn = model.recede(&unloaded.candidate, &recession, cx).unwrap();
            let mut softened = config.clone();
            softened.hardness_pa[12] *= 0.5;
            let successor =
                FiniteGapPlasticModel::new_with_linear_solver(softened, selection, cx).unwrap();
            let refreshed = successor.refresh_from(&model, &worn.candidate, cx).unwrap();
            assert_eq!(
                refreshed.candidate.plastic_depth_m(),
                step.candidate.plastic_depth_m()
            );
            assert_eq!(refreshed.candidate.recession_m(), recession);
            assert_eq!(
                refreshed.candidate.cumulative_plastic_work_j(),
                step.candidate.cumulative_plastic_work_j()
            );
        }
    });
}

#[test]
fn iterative_work_limits_ownership_and_cancellation_refuse_without_mutation() {
    with_cx(|cx, gate| {
        let mut config = small_config();
        config.undeformed_gap_m[11] = 0.3e-6;
        config.hardness_pa.fill(1.0e9);
        let limited = FiniteGapPlasticModel::new_with_linear_solver(
            config.clone(),
            FiniteGapPlasticLinearSolver::MatrixFree {
                maximum_iterations: 1,
                maximum_corrections: 0,
            },
            cx,
        )
        .unwrap();
        let before = limited.initial_state();
        let saved = before.clone();
        assert!(
            matches!(limited.advance(&before, 2.0e-6, cx), Err(FiniteGapPlasticError::LinearSolveDidNotConverge {
            free_cells: 3, iterations: 1, residual_m,
        }) if residual_m.is_finite() && residual_m > config.complementarity_tolerance_m)
        );
        assert_eq!(before, saved);
        let full =
            FiniteGapPlasticModel::new_with_linear_solver(config.clone(), ITERATIVE, cx).unwrap();
        let dense = FiniteGapPlasticModel::new(config, cx).unwrap();
        assert_ne!(full.identity(), dense.identity());
        assert_eq!(
            full.advance(&before, 2.0e-6, cx),
            Err(FiniteGapPlasticError::ForeignState)
        );
        assert!(matches!(
            full.refresh_from(&dense, &dense.initial_state(), cx),
            Err(FiniteGapPlasticError::Contact(
                FiniteGapContactError::InvalidInput { .. }
            ))
        ));
        let before = full.initial_state();
        let saved = before.clone();
        assert!(full.advance(&before, 2.0e-6, cx).is_ok());
        gate.request();
        assert_eq!(
            full.advance(&before, 2.0e-6, cx),
            Err(FiniteGapPlasticError::Contact(
                FiniteGapContactError::Cancelled
            ))
        );
        assert_eq!(before, saved);
    });
}
