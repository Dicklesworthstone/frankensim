//! G0/G1/G3 local derivatives of the hardness-limited return map.
//! Re-solves are finite-difference oracles, not experimental material validation.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::normal_patch::{
    FiniteGapGrid, FiniteGapPlasticConfig, FiniteGapPlasticDirection, FiniteGapPlasticError,
    FiniteGapPlasticLinearSolver, FiniteGapPlasticModel, FiniteGapPlasticSensitivityControls,
    FiniteGapPlasticSensitivityError,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx<T>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 81,
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
        maximum_active_set_iterations: 512,
        complementarity_tolerance_m: 1e-14,
        boundary_clearance_cells: 1,
    }
}

fn controls() -> FiniteGapPlasticSensitivityControls {
    FiniteGapPlasticSensitivityControls {
        minimum_pressure_margin_pa: 100.0,
        minimum_displacement_margin_m: 1e-10,
        maximum_linear_residual_m: 1e-13,
    }
}

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual {actual:.16e}, expected {expected:.16e}"
    );
}

fn self_compliance() -> f64 {
    // Independent square-cell integral of 1/(pi E* r).
    4.0 * 1e-3 * (1.0 + 2.0_f64.sqrt()).ln() / (core::f64::consts::PI * 1e9)
}

#[test]
fn independent_elastic_stiffness_and_plastic_plateau_tangents() {
    with_cx(|cx, _| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let before = model.initial_state();
        let zero = vec![0.0; 25];
        let c = self_compliance();
        let h = model.config().hardness_pa[12];
        let area = 1e-6;
        let direction = FiniteGapPlasticDirection {
            approach_m: 1.0,
            reduced_modulus_pa: 0.0,
            hardness_pa: &zero,
            recession_m: &zero,
        };
        let elastic = model.advance(&before, 0.5 * c * h, cx).unwrap();
        let tangent = model
            .linearize_advance(&before, &elastic, controls(), cx)
            .unwrap();
        assert_eq!(tangent.margins().elastic_cells, 1);
        assert_eq!(tangent.margins().flowing_cells, 0);
        let derivative = tangent.apply(direction, cx).unwrap();
        close(derivative.pressure_pa[12], 1.0 / c, 2e-13, 0.0);
        close(derivative.normal_force_n, area / c, 2e-13, 0.0);
        close(
            derivative.elastic_energy_j,
            elastic.receipt.normal.normal_force_n,
            2e-13,
            0.0,
        );
        assert_eq!(derivative.plastic_depth_m, zero);
        assert_eq!(derivative.plastic_work_j, 0.0);

        let flowing = model.advance(&before, c * h + 3e-6, cx).unwrap();
        let tangent = model
            .linearize_advance(&before, &flowing, controls(), cx)
            .unwrap();
        assert_eq!(tangent.margins().flowing_cells, 1);
        assert_eq!(tangent.margins().elastic_cells, 0);
        let derivative = tangent.apply(direction, cx).unwrap();
        assert_eq!(derivative.pressure_pa, zero);
        assert_eq!(derivative.normal_force_n, 0.0);
        close(derivative.plastic_depth_m[12], 1.0, 2e-14, 0.0);
        close(derivative.plastic_work_j, h * area, 2e-14, 0.0);
        assert_eq!(derivative.elastic_energy_j, 0.0);
        assert_eq!(derivative.model_identity, model.identity());
        assert_eq!(derivative.predecessor_identity, before.identity());
        assert_eq!(derivative.candidate_identity, flowing.candidate.identity());

        let mut dh = zero.clone();
        dh[12] = -1.0;
        let hardness = tangent
            .apply(
                FiniteGapPlasticDirection {
                    approach_m: 0.0,
                    hardness_pa: &dh,
                    ..direction
                },
                cx,
            )
            .unwrap();
        close(hardness.pressure_pa[12], -1.0, 0.0, 0.0);
        close(hardness.elastic_opening_m[12], -c, 2e-14, 0.0);
        close(hardness.plastic_depth_m[12], c, 2e-14, 0.0);
        close(hardness.normal_force_n, -area, 2e-14, 0.0);
        close(
            hardness.plastic_work_j,
            (h * c - 3e-6) * area,
            2e-13,
            0.0,
        );

        let modulus = tangent
            .apply(
                FiniteGapPlasticDirection {
                    approach_m: 0.0,
                    reduced_modulus_pa: model.config().reduced_modulus_pa,
                    ..direction
                },
                cx,
            )
            .unwrap();
        assert_eq!(modulus.pressure_pa, zero);
        close(modulus.elastic_opening_m[12], -c * h, 2e-14, 0.0);
        close(modulus.plastic_depth_m[12], c * h, 2e-14, 0.0);
        close(
            modulus.elastic_energy_j,
            -flowing.receipt.normal.reversible_energy_j,
            2e-13,
            0.0,
        );
        close(modulus.plastic_work_j, h * area * c * h, 2e-13, 0.0);
    });
}

#[test]
fn mixed_active_sets_match_actual_material_geometry_and_approach_resolves() {
    with_cx(|cx, _| {
        for solver in [
            FiniteGapPlasticLinearSolver::Dense,
            FiniteGapPlasticLinearSolver::MatrixFree {
                maximum_iterations: 256,
                maximum_corrections: 3,
            },
            FiniteGapPlasticLinearSolver::Fft {
                maximum_iterations: 256,
                maximum_corrections: 3,
            },
        ] {
            let mut cfg = config();
            cfg.undeformed_gap_m[12] = 0.01;
            cfg.undeformed_gap_m[11] = 0.0;
            cfg.undeformed_gap_m[13] = 0.0;
            cfg.hardness_pa[13] = 3e6;
            let model =
                FiniteGapPlasticModel::new_with_linear_solver(cfg.clone(), solver, cx).unwrap();
            let mut recession = vec![0.0; 25];
            recession[11] = 1e-6;
            recession[13] = 1e-6;
            let before = model
                .recede(&model.initial_state(), &recession, cx)
                .unwrap()
                .candidate;
            let approach = 3e-6;
            let advance = model.advance(&before, approach, cx).unwrap();
            let tangent = model
                .linearize_advance(&before, &advance, controls(), cx)
                .unwrap();
            assert_eq!(tangent.margins().elastic_cells, 1);
            assert_eq!(tangent.margins().flowing_cells, 1);
            assert_eq!(tangent.margins().open_cells, 23);
            let mut dh = vec![0.0; 25];
            dh[11] = -2e5;
            dh[13] = 3e5;
            let mut dr = vec![0.0; 25];
            dr[11] = 2e-7;
            dr[13] = -1e-7;
            let da = 4e-7;
            let de = 2e8;
            let derivative = tangent
                .apply(
                    FiniteGapPlasticDirection {
                        approach_m: da,
                        reduced_modulus_pa: de,
                        hardness_pa: &dh,
                        recession_m: &dr,
                    },
                    cx,
                )
                .unwrap();
            assert!(derivative.pressure_pa[11] < 0.0);
            assert!(derivative.linear_residual_m <= controls().maximum_linear_residual_m);

            // Each endpoint solves the actual native law from the same
            // physical zero-plastic-history state, with perturbed recession.
            let resolve = |epsilon: f64| {
                let mut perturbed = cfg.clone();
                perturbed.reduced_modulus_pa += epsilon * de;
                for (h, delta) in perturbed.hardness_pa.iter_mut().zip(&dh) {
                    *h += epsilon * delta;
                }
                let owner =
                    FiniteGapPlasticModel::new_with_linear_solver(perturbed, solver, cx).unwrap();
                let profile: Vec<f64> = recession
                    .iter()
                    .zip(&dr)
                    .map(|(r, delta)| r + epsilon * delta)
                    .collect();
                let retained = owner
                    .recede(&owner.initial_state(), &profile, cx)
                    .unwrap()
                    .candidate;
                owner
                    .advance(&retained, approach + epsilon * da, cx)
                    .unwrap()
            };
            let epsilon = 1e-3;
            let plus = resolve(epsilon);
            let minus = resolve(-epsilon);
            let difference = |p: f64, m: f64| (p - m) / (2.0 * epsilon);
            for i in 0..25 {
                close(
                    derivative.pressure_pa[i],
                    difference(
                        plus.candidate.pressure_pa()[i],
                        minus.candidate.pressure_pa()[i],
                    ),
                    2e-6,
                    1e-3,
                );
                close(
                    derivative.plastic_depth_m[i],
                    difference(
                        plus.candidate.plastic_depth_m()[i],
                        minus.candidate.plastic_depth_m()[i],
                    ),
                    2e-6,
                    1e-13,
                );
                close(
                    derivative.elastic_opening_m[i],
                    difference(
                        plus.receipt.normal.elastic_opening_m[i],
                        minus.receipt.normal.elastic_opening_m[i],
                    ),
                    2e-6,
                    1e-13,
                );
            }
            close(
                derivative.normal_force_n,
                difference(
                    plus.receipt.normal.normal_force_n,
                    minus.receipt.normal.normal_force_n,
                ),
                2e-6,
                1e-9,
            );
            close(
                derivative.elastic_energy_j,
                difference(
                    plus.receipt.normal.reversible_energy_j,
                    minus.receipt.normal.reversible_energy_j,
                ),
                2e-6,
                1e-13,
            );
            close(
                derivative.plastic_work_j,
                difference(plus.receipt.plastic_work_j, minus.receipt.plastic_work_j),
                2e-6,
                1e-13,
            );
        }
    });
}

#[test]
fn opening_yield_and_unloading_corners_refuse_a_classical_tangent() {
    with_cx(|cx, _| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let before = model.initial_state();
        let touch = model.advance(&before, 0.0, cx).unwrap();
        assert!(matches!(
            model.linearize_advance(&before, &touch, controls(), cx),
            Err(FiniteGapPlasticSensitivityError::SwitchingCell { cell: 12, .. })
        ));
        let yield_approach = self_compliance() * 1e6;
        let yielding = model.advance(&before, yield_approach, cx).unwrap();
        assert!(matches!(
            model.linearize_advance(&before, &yielding, controls(), cx),
            Err(FiniteGapPlasticSensitivityError::SwitchingCell { cell: 12, .. })
        ));
        let flow = model
            .advance(&before, yield_approach + 3e-6, cx)
            .unwrap();
        let corner = model
            .advance(&flow.candidate, flow.candidate.approach_m(), cx)
            .unwrap();
        assert!(matches!(
            model.linearize_advance(&flow.candidate, &corner, controls(), cx),
            Err(FiniteGapPlasticSensitivityError::SwitchingCell { cell: 12, .. })
        ));

        let unload = model
            .advance(
                &flow.candidate,
                flow.candidate.approach_m() - 0.25 * yield_approach,
                cx,
            )
            .unwrap();
        let tangent = model
            .linearize_advance(&flow.candidate, &unload, controls(), cx)
            .unwrap();
        let zero = vec![0.0; 25];
        let derivative = tangent
            .apply(
                FiniteGapPlasticDirection {
                    approach_m: 1.0,
                    reduced_modulus_pa: 0.0,
                    hardness_pa: &zero,
                    recession_m: &zero,
                },
                cx,
            )
            .unwrap();
        close(
            derivative.normal_force_n,
            1e-6 / self_compliance(),
            2e-13,
            0.0,
        );
        assert_eq!(derivative.plastic_depth_m, zero);
        assert_eq!(derivative.plastic_work_j, 0.0);
        assert_eq!(
            unload.candidate.plastic_depth_m(),
            flow.candidate.plastic_depth_m(),
        );
    });
}

#[test]
fn ownership_direction_controls_and_cancellation_leave_history_untouched() {
    with_cx(|cx, gate| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let before = model.initial_state();
        let advance = model
            .advance(&before, 0.5 * self_compliance() * 1e6, cx)
            .unwrap();
        let saved = advance.clone();
        let mut other_cfg = config();
        other_cfg.hardness_pa[12] *= 2.0;
        let other = FiniteGapPlasticModel::new(other_cfg, cx).unwrap();
        assert!(matches!(
            other.linearize_advance(&before, &advance, controls(), cx),
            Err(FiniteGapPlasticSensitivityError::Contact(
                FiniteGapPlasticError::ForeignState
            ))
        ));
        let mut invalid = controls();
        invalid.minimum_displacement_margin_m = model.config().complementarity_tolerance_m;
        assert!(
            model
                .linearize_advance(&before, &advance, invalid, cx)
                .is_err()
        );
        let tangent = model
            .linearize_advance(&before, &advance, controls(), cx)
            .unwrap();
        let zero = vec![0.0; 25];
        let direction = FiniteGapPlasticDirection {
            approach_m: 0.0,
            reduced_modulus_pa: 0.0,
            hardness_pa: &zero,
            recession_m: &zero,
        };
        let result = tangent.apply(direction, cx).unwrap();
        assert_eq!(result.pressure_pa, zero);
        assert_eq!(result.elastic_energy_j, 0.0);
        assert!(
            tangent
                .apply(
                    FiniteGapPlasticDirection {
                        hardness_pa: &zero[..24],
                        ..direction
                    },
                    cx,
                )
                .is_err()
        );
        assert!(
            tangent
                .apply(
                    FiniteGapPlasticDirection {
                        reduced_modulus_pa: f64::NAN,
                        ..direction
                    },
                    cx,
                )
                .is_err()
        );
        let mut wrong_predecessor = advance.clone();
        wrong_predecessor.receipt.predecessor_identity = other.initial_state().identity();
        assert!(
            model
                .linearize_advance(&before, &wrong_predecessor, controls(), cx)
                .is_err()
        );
        gate.request();
        assert!(tangent.apply(direction, cx).is_err());
        assert_eq!(advance, saved);
    });
}
