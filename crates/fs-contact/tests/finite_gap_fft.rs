//! Physical references for isolated-domain FFT pressure solves. These are
//! numerical contact checks, not material calibration or timing measurements.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::normal_patch::{
    FiniteGapContactError, FiniteGapGrid, FiniteGapPlasticConfig, FiniteGapPlasticDrive,
    FiniteGapPlasticDriveAdvance, FiniteGapPlasticError, FiniteGapPlasticLinearSolver,
    FiniteGapPlasticLoadControls, FiniteGapPlasticModel, FiniteGapPlasticReceipt,
    FiniteGapPlasticState, MAX_FINITE_GAP_PLASTIC_FREE_CELLS,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

const FFT: FiniteGapPlasticLinearSolver = FiniteGapPlasticLinearSolver::Fft {
    maximum_iterations: 512,
    maximum_corrections: 2,
};
const SELECTIONS: [FiniteGapPlasticLinearSolver; 3] = [
    FiniteGapPlasticLinearSolver::Dense,
    FiniteGapPlasticLinearSolver::MatrixFree {
        maximum_iterations: 512,
        maximum_corrections: 2,
    },
    FFT,
];

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 193,
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

// Independent rectangular 1/r integral, using asinh rather than the native
// logarithmic primitive. Both physical cell dimensions remain explicit.
fn rectangle_compliance(x: f64, y: f64, grid: FiniteGapGrid, modulus: f64) -> f64 {
    fn primitive(x: f64, y: f64) -> f64 {
        let first = if x == 0.0 {
            0.0
        } else {
            x * (y / x.abs()).asinh()
        };
        let second = if y == 0.0 {
            0.0
        } else {
            y * (x / y.abs()).asinh()
        };
        first + second
    }
    let hx = 0.5 * grid.cell_width_m;
    let hy = 0.5 * grid.cell_depth_m;
    (primitive(x + hx, y + hy) - primitive(x - hx, y + hy) - primitive(x + hx, y - hy)
        + primitive(x - hx, y - hy))
        / (core::f64::consts::PI * modulus)
}

struct MixedPatch {
    config: FiniteGapPlasticConfig,
    approach_m: f64,
    pressure_pa: Vec<f64>,
    opening_m: Vec<f64>,
    plastic_depth_m: Vec<f64>,
}

fn mixed_patch() -> MixedPatch {
    let grid = FiniteGapGrid {
        cells_x: 7,
        cells_y: 5,
        cell_width_m: 0.8e-3,
        cell_depth_m: 1.3e-3,
    };
    let count = grid.cells_x * grid.cells_y;
    let modulus = 1e9;
    let mut pressure = vec![0.0; count];
    for (index, value) in [(9, 0.25e6), (18, 0.9e6), (22, 0.5e6), (12, 0.6e6)] {
        pressure[index] = value;
    }
    let mut hardness = vec![1.8e6; count];
    hardness[9] = pressure[9];
    hardness[22] = pressure[22];
    let mut plastic_depth = vec![0.0; count];
    plastic_depth[9] = 0.4e-6;
    plastic_depth[22] = 0.15e-6;
    let mut opening = vec![0.0; count];
    let mut gap = vec![0.01; count];
    let mut approach = 0.0_f64;
    for i in 0..count {
        for (j, &p) in pressure.iter().enumerate().filter(|(_, p)| **p > 0.0) {
            let dx = (i % grid.cells_x).abs_diff(j % grid.cells_x) as f64;
            let dy = (i / grid.cells_x).abs_diff(j / grid.cells_x) as f64;
            opening[i] += rectangle_compliance(
                dx * grid.cell_width_m,
                dy * grid.cell_depth_m,
                grid,
                modulus,
            ) * p;
        }
        if pressure[i] > 0.0 {
            gap[i] = -opening[i] - plastic_depth[i];
            approach = approach.max(-gap[i]);
        }
    }
    MixedPatch {
        config: FiniteGapPlasticConfig {
            grid,
            undeformed_gap_m: gap,
            reduced_modulus_pa: modulus,
            hardness_pa: hardness,
            maximum_active_set_iterations: 128,
            complementarity_tolerance_m: 1e-13,
            boundary_clearance_cells: 1,
        },
        approach_m: approach,
        pressure_pa: pressure,
        opening_m: opening,
        plastic_depth_m: plastic_depth,
    }
}

fn assert_same_physics(
    actual: &FiniteGapPlasticState,
    response: &FiniteGapPlasticReceipt,
    expected: &FiniteGapPlasticState,
    reference: &FiniteGapPlasticReceipt,
) {
    close(actual.approach_m(), expected.approach_m(), 0.0, 1e-13);
    for (values, originals, relative, absolute) in [
        (actual.pressure_pa(), expected.pressure_pa(), 2e-7, 1e-6),
        (
            actual.plastic_depth_m(),
            expected.plastic_depth_m(),
            0.0,
            1e-13,
        ),
        (actual.recession_m(), expected.recession_m(), 0.0, 1e-15),
        (
            response.normal.elastic_opening_m.as_slice(),
            reference.normal.elastic_opening_m.as_slice(),
            0.0,
            1e-13,
        ),
    ] {
        assert_eq!(values.len(), originals.len());
        for (&value, &original) in values.iter().zip(originals) {
            close(value, original, relative, absolute);
        }
    }
    assert_eq!(response.normal.active_cells, reference.normal.active_cells);
    close(
        response.normal.normal_force_n,
        reference.normal.normal_force_n,
        2e-7,
        1e-12,
    );
    for (value, original) in [
        (
            response.normal.reversible_energy_j,
            reference.normal.reversible_energy_j,
        ),
        (response.plastic_work_j, reference.plastic_work_j),
        (
            actual.cumulative_plastic_work_j(),
            expected.cumulative_plastic_work_j(),
        ),
    ] {
        close(value, original, 2e-7, 1e-18);
    }
    assert!(response.kkt_residual_m <= 1e-13);
    assert!(response.normal.complementarity_residual_m <= 1e-13);
}

#[test]
fn rectangular_off_axis_mixed_caps_preserve_physics_through_contact_history() {
    with_cx(|cx, _| {
        let fixture = mixed_patch();
        let mut reference: Option<Vec<(FiniteGapPlasticState, FiniteGapPlasticReceipt)>> = None;
        for selection in SELECTIONS {
            let model = FiniteGapPlasticModel::new_with_linear_solver(
                fixture.config.clone(),
                selection,
                cx,
            )
            .unwrap();
            let initial = model.initial_state();
            let loaded = model.advance(&initial, fixture.approach_m, cx).unwrap();
            assert_eq!(loaded.receipt.yielded_cells, 2);
            assert_eq!(loaded.receipt.normal.active_cells, 4);
            for i in 0..fixture.pressure_pa.len() {
                close(
                    loaded.candidate.pressure_pa()[i],
                    fixture.pressure_pa[i],
                    2e-7,
                    1e-6,
                );
                close(
                    loaded.candidate.plastic_depth_m()[i],
                    fixture.plastic_depth_m[i],
                    0.0,
                    1e-13,
                );
                close(
                    loaded.receipt.normal.elastic_opening_m[i],
                    fixture.opening_m[i],
                    0.0,
                    1e-13,
                );
            }
            assert_eq!(initial, model.initial_state());
            let unloaded = model.advance(&loaded.candidate, -1e-6, cx).unwrap();
            assert_eq!(unloaded.receipt.normal.normal_force_n, 0.0);
            assert_eq!(
                unloaded.candidate.plastic_depth_m(),
                loaded.candidate.plastic_depth_m()
            );
            assert_eq!(
                unloaded.candidate.cumulative_plastic_work_j(),
                loaded.candidate.cumulative_plastic_work_j()
            );
            let reloaded = model
                .advance(&unloaded.candidate, fixture.approach_m, cx)
                .unwrap();
            assert!(reloaded.receipt.plastic_work_j <= 1e-15);
            let mut recession = vec![0.0; fixture.pressure_pa.len()];
            recession[9] = 0.15e-6;
            recession[18] = 0.1e-6;
            let worn = model.recede(&reloaded.candidate, &recession, cx).unwrap();
            assert_eq!(worn.candidate.recession_m(), recession);
            assert_eq!(worn.receipt.normal_work_j, 0.0);
            close(
                worn.receipt.parameter_work_j,
                worn.receipt.elastic_energy_change_j + worn.receipt.response.plastic_work_j,
                1e-12,
                1e-20,
            );
            let saved_worn = worn.candidate.clone();
            let mut changed = fixture.config.clone();
            changed.reduced_modulus_pa *= 0.8;
            changed.hardness_pa[18] = 0.3e6;
            changed.hardness_pa[9] *= 0.6;
            let successor =
                FiniteGapPlasticModel::new_with_linear_solver(changed, selection, cx).unwrap();
            let refreshed = successor.refresh_from(&model, &worn.candidate, cx).unwrap();
            assert_eq!(successor.linear_solver(), selection);
            assert_eq!(refreshed.candidate.owner_identity(), successor.identity());
            assert_eq!(refreshed.candidate.recession_m(), recession);
            assert_eq!(
                refreshed.receipt.response.predecessor_identity,
                worn.candidate.identity()
            );
            assert!(refreshed.receipt.response.plastic_work_j > 0.0);
            assert_eq!(refreshed.receipt.normal_work_j, 0.0);
            close(
                refreshed.receipt.parameter_work_j,
                refreshed.receipt.elastic_energy_change_j
                    + refreshed.receipt.response.plastic_work_j,
                1e-12,
                1e-20,
            );
            assert_eq!(worn.candidate, saved_worn);
            let snapshots = vec![
                (loaded.candidate, loaded.receipt),
                (unloaded.candidate, unloaded.receipt),
                (reloaded.candidate, reloaded.receipt),
                (worn.candidate, worn.receipt.response),
                (refreshed.candidate, refreshed.receipt.response),
            ];
            if let Some(reference) = &reference {
                for ((state, response), (previous, receipt)) in snapshots.iter().zip(reference) {
                    assert_same_physics(state, response, previous, receipt);
                }
            } else {
                reference = Some(snapshots);
            }
        }
    });
}

#[test]
fn rank_one_spring_coupling_retains_native_pressure_work_and_contact_storage() {
    with_cx(|cx, _| {
        let fixture = mixed_patch();
        let area = fixture.config.grid.cell_width_m * fixture.config.grid.cell_depth_m;
        let force = fixture.pressure_pa.iter().sum::<f64>() * area;
        let stiffness = 4e5;
        let drive = FiniteGapPlasticDrive {
            stiffness_n_per_m: stiffness,
            position_m: fixture.approach_m + force / stiffness,
        };
        let controls = FiniteGapPlasticLoadControls {
            minimum_approach_m: -0.01,
            maximum_approach_m: 0.02,
            absolute_force_tolerance_n: 1e-12,
            relative_force_tolerance: 1e-12,
            approach_tolerance_m: 1e-18,
            maximum_equilibrium_evaluations: 1,
        };
        let elastic_energy = 0.5
            * area
            * fixture
                .pressure_pa
                .iter()
                .zip(&fixture.opening_m)
                .map(|(p, opening)| p * opening)
                .sum::<f64>();
        let plastic_work = area
            * fixture
                .pressure_pa
                .iter()
                .zip(&fixture.plastic_depth_m)
                .map(|(p, dent)| p * dent)
                .sum::<f64>();
        let mut reference: Option<FiniteGapPlasticDriveAdvance> = None;
        for selection in SELECTIONS {
            let model = FiniteGapPlasticModel::new_with_linear_solver(
                fixture.config.clone(),
                selection,
                cx,
            )
            .unwrap();
            let initial = model.initial_state();
            let next = model
                .advance_with_drive(&initial, drive, controls, cx)
                .unwrap();
            close(next.drive.approach_m, fixture.approach_m, 0.0, 1e-13);
            close(next.drive.achieved_force_n, force, 2e-7, 1e-12);
            close(
                next.contact.receipt.normal.reversible_energy_j,
                elastic_energy,
                2e-7,
                1e-18,
            );
            close(
                next.contact.receipt.plastic_work_j,
                plastic_work,
                2e-7,
                1e-18,
            );
            close(
                next.drive.spring_energy_j,
                0.5 * force * force / stiffness,
                2e-7,
                1e-18,
            );
            assert_eq!(
                next.drive.achieved_force_n,
                next.contact.candidate.pressure_pa().iter().sum::<f64>() * area
            );
            assert_eq!(
                next.drive.drive_force_n,
                stiffness * (drive.position_m - next.drive.approach_m)
            );
            assert!(next.drive.force_residual_n.abs() <= next.drive.force_tolerance_n);
            assert_eq!(next.drive.equilibrium_evaluations, 1);
            assert_eq!(next.drive.approach_bracket_m, [next.drive.approach_m; 2]);
            assert_eq!(initial, model.initial_state());
            if let Some(reference) = &reference {
                assert_same_physics(
                    &next.contact.candidate,
                    &next.contact.receipt,
                    &reference.contact.candidate,
                    &reference.contact.receipt,
                );
            } else {
                reference = Some(next);
            }
        }
    });
}

#[test]
fn fft_iteration_limits_owner_selection_and_cancellation_remain_binding() {
    with_cx(|cx, gate| {
        let fixture = mixed_patch();
        let mut config = fixture.config;
        config.hardness_pa.fill(1e9);
        let limited = FiniteGapPlasticModel::new_with_linear_solver(
            config.clone(),
            FiniteGapPlasticLinearSolver::Fft {
                maximum_iterations: 1,
                maximum_corrections: 0,
            },
            cx,
        )
        .unwrap();
        let before = limited.initial_state();
        let saved = before.clone();
        assert!(matches!(limited.advance(&before, fixture.approach_m, cx),
            Err(FiniteGapPlasticError::LinearSolveDidNotConverge {
                free_cells: 4, iterations: 1, residual_m,
            }) if residual_m.is_finite() && residual_m > config.complementarity_tolerance_m));
        assert_eq!(before, saved);
        assert!(matches!(
            FiniteGapPlasticModel::new_with_linear_solver(
                config.clone(),
                FiniteGapPlasticLinearSolver::Fft {
                    maximum_iterations: 0,
                    maximum_corrections: 2
                },
                cx,
            ),
            Err(FiniteGapPlasticError::Contact(
                FiniteGapContactError::InvalidInput { .. }
            ))
        ));
        let model = FiniteGapPlasticModel::new_with_linear_solver(config.clone(), FFT, cx).unwrap();
        let dense = FiniteGapPlasticModel::new(config.clone(), cx).unwrap();
        let direct_cg =
            FiniteGapPlasticModel::new_with_linear_solver(config.clone(), SELECTIONS[1], cx)
                .unwrap();
        assert_ne!(model.identity(), dense.identity());
        assert_ne!(model.identity(), direct_cg.identity());
        assert_ne!(model.identity(), limited.identity());
        assert_eq!(
            model.advance(&before, fixture.approach_m, cx),
            Err(FiniteGapPlasticError::ForeignState)
        );
        assert!(matches!(
            model.refresh_from(&dense, &dense.initial_state(), cx),
            Err(FiniteGapPlasticError::Contact(
                FiniteGapContactError::InvalidInput { .. }
            ))
        ));
        let before = model.initial_state();
        let saved = before.clone();
        assert!(model.advance(&before, fixture.approach_m, cx).is_ok());
        gate.request();
        assert_eq!(
            model.advance(&before, fixture.approach_m, cx),
            Err(FiniteGapPlasticError::Contact(
                FiniteGapContactError::Cancelled
            ))
        );
        assert!(matches!(
            FiniteGapPlasticModel::new_with_linear_solver(config, FFT, cx),
            Err(FiniteGapPlasticError::Contact(
                FiniteGapContactError::Cancelled
            ))
        ));
        assert_eq!(before, saved);
    });
}

#[test]
fn manufactured_rectangular_patch_exceeds_dense_limit_without_periodic_contact() {
    with_cx(|cx, _| {
        let grid = FiniteGapGrid {
            cells_x: 49,
            cells_y: 47,
            cell_width_m: 1.1e-3,
            cell_depth_m: 0.7e-3,
        };
        let count = grid.cells_x * grid.cells_y;
        let modulus = 1e9;
        let active: Vec<usize> = (0..count)
            .filter(|i| {
                let x = i % grid.cells_x;
                let y = i / grid.cells_x;
                x > 0 && x < grid.cells_x - 1 && y > 0 && y < grid.cells_y - 1
            })
            .collect();
        assert_eq!(active.len(), 2115);
        assert!(active.len() > MAX_FINITE_GAP_PLASTIC_FREE_CELLS);
        let kernel: Vec<f64> = (0..count)
            .map(|i| {
                rectangle_compliance(
                    (i % grid.cells_x) as f64 * grid.cell_width_m,
                    (i / grid.cells_x) as f64 * grid.cell_depth_m,
                    grid,
                    modulus,
                )
            })
            .collect();
        let mut pressure = vec![0.0; count];
        for &i in &active {
            let x = i % grid.cells_x;
            let y = i / grid.cells_x;
            pressure[i] =
                8e4 + 370.0 * x as f64 + 910.0 * y as f64 + 2500.0 * ((x + 2 * y) % 3) as f64;
        }
        let mut opening = vec![0.0; count];
        let mut gap = vec![1.0; count];
        let mut approach = 0.0_f64;
        for &i in &active {
            opening[i] = active
                .iter()
                .map(|&j| {
                    let dx = (i % grid.cells_x).abs_diff(j % grid.cells_x);
                    let dy = (i / grid.cells_x).abs_diff(j / grid.cells_x);
                    kernel[dy * grid.cells_x + dx] * pressure[j]
                })
                .sum();
            gap[i] = -opening[i];
            approach = approach.max(opening[i]);
        }
        let config = FiniteGapPlasticConfig {
            grid,
            undeformed_gap_m: gap,
            reduced_modulus_pa: modulus,
            hardness_pa: vec![1e9; count],
            maximum_active_set_iterations: 32,
            complementarity_tolerance_m: 1e-11,
            boundary_clearance_cells: 1,
        };
        let dense = FiniteGapPlasticModel::new(config.clone(), cx).unwrap();
        assert!(matches!(
            dense.advance(&dense.initial_state(), approach, cx),
            Err(FiniteGapPlasticError::FreeSetTooLarge {
                free_cells: 2115,
                maximum: 2048
            })
        ));
        let model = FiniteGapPlasticModel::new_with_linear_solver(config.clone(), FFT, cx).unwrap();
        let before = model.initial_state();
        let saved = before.clone();
        let next = model.advance(&before, approach, cx).unwrap();
        let area = grid.cell_width_m * grid.cell_depth_m;
        let expected_force = pressure.iter().sum::<f64>() * area;
        let expected_energy = 0.5
            * area
            * pressure
                .iter()
                .zip(&opening)
                .map(|(p, u)| p * u)
                .sum::<f64>();
        assert_eq!(next.receipt.normal.active_cells, active.len());
        assert_eq!(next.receipt.yielded_cells, 0);
        assert_eq!(next.receipt.plastic_work_j, 0.0);
        assert_eq!(next.candidate.cumulative_plastic_work_j(), 0.0);
        assert!(next.candidate.plastic_depth_m().iter().all(|z| *z == 0.0));
        assert!(next.receipt.kkt_residual_m <= config.complementarity_tolerance_m);
        assert!(
            next.receipt.normal.complementarity_residual_m <= config.complementarity_tolerance_m
        );
        close(
            next.receipt.normal.normal_force_n,
            expected_force,
            5e-6,
            0.0,
        );
        close(
            next.receipt.normal.reversible_energy_j,
            expected_energy,
            5e-6,
            0.0,
        );
        for (i, &expected) in pressure.iter().enumerate() {
            if expected == 0.0 {
                assert_eq!(next.candidate.pressure_pa()[i], 0.0);
            } else {
                close(next.candidate.pressure_pa()[i], expected, 1e-4, 0.0);
                close(
                    next.receipt.normal.elastic_opening_m[i],
                    opening[i],
                    0.0,
                    config.complementarity_tolerance_m,
                );
            }
        }
        assert_eq!(before, saved);
    });
}
