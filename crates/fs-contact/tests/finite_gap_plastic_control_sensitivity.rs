//! G0/G1/G3 fixed-history tangents under force and finite-spring boundaries.
//! Independent scalar references and actual re-solves are numerical-model checks.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_contact::normal_patch::{
    FiniteGapGrid, FiniteGapPlasticAdvance, FiniteGapPlasticConfig, FiniteGapPlasticDrive,
    FiniteGapPlasticDriveDirection, FiniteGapPlasticDriveSensitivityControls,
    FiniteGapPlasticLinearSolver, FiniteGapPlasticLoadControls, FiniteGapPlasticLoadDirection,
    FiniteGapPlasticLoadError, FiniteGapPlasticLoadSensitivityControls, FiniteGapPlasticModel,
    FiniteGapPlasticSensitivity, FiniteGapPlasticSensitivityControls,
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
                seed: 73,
                kernel_id: 83,
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

fn equilibrium_controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1e-5,
        maximum_approach_m: 2e-5,
        absolute_force_tolerance_n: 1e-13,
        relative_force_tolerance: 1e-12,
        approach_tolerance_m: 1e-18,
        maximum_equilibrium_evaluations: 128,
    }
}

fn contact_controls() -> FiniteGapPlasticSensitivityControls {
    FiniteGapPlasticSensitivityControls {
        minimum_pressure_margin_pa: 100.0,
        minimum_displacement_margin_m: 1e-10,
        maximum_linear_residual_m: 1e-13,
    }
}

fn load_controls() -> FiniteGapPlasticLoadSensitivityControls {
    FiniteGapPlasticLoadSensitivityControls {
        contact: contact_controls(),
        minimum_contact_stiffness_n_per_m: 1.0,
        maximum_force_residual_n: 1e-10,
    }
}

fn drive_controls() -> FiniteGapPlasticDriveSensitivityControls {
    FiniteGapPlasticDriveSensitivityControls {
        contact: contact_controls(),
        maximum_force_residual_n: 1e-10,
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
fn force_tangent_retains_the_dent_and_matches_independent_elastic_balance() {
    with_cx(|cx, _| {
        let c = self_compliance();
        let area = 1e-6;
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let before = model
            .advance(&model.initial_state(), c * 1e6 + 2e-6, cx)
            .unwrap()
            .candidate;
        let saved = before.clone();
        let advance = model
            .advance_at_load(&before, 0.4, equilibrium_controls(), cx)
            .unwrap();
        let tangent = model
            .linearize_load_advance(&before, &advance, load_controls(), cx)
            .unwrap();
        close(tangent.contact_stiffness_n_per_m(), area / c, 2e-13, 0.0);
        assert_eq!(tangent.margins().elastic_cells, 1);
        let dh = vec![-3e5; 25];
        let mut dr = vec![0.0; 25];
        dr[12] = 2e-7;
        let direction = FiniteGapPlasticLoadDirection {
            target_force_n: 0.2,
            reduced_modulus_pa: 2e8,
            hardness_pa: &dh,
            recession_m: &dr,
        };
        let derivative = tangent.apply(direction, cx).unwrap();
        let p = 0.4 / area;
        let dp = direction.target_force_n / area;
        close(derivative.contact.pressure_pa[12], dp, 2e-12, 1e-8);
        close(
            derivative.approach_m,
            dr[12] + c * (dp - 0.2 * p),
            2e-12,
            1e-18,
        );
        close(derivative.contact.normal_force_n, 0.2, 2e-12, 1e-14);
        close(
            derivative.contact.elastic_energy_j,
            area * c * p * (dp - 0.1 * p),
            2e-12,
            1e-18,
        );
        assert_eq!(derivative.contact.plastic_depth_m, vec![0.0; 25]);
        assert_eq!(derivative.contact.plastic_work_j, 0.0);
        let maintained = tangent
            .apply(
                FiniteGapPlasticLoadDirection {
                    target_force_n: 0.0,
                    ..direction
                },
                cx,
            )
            .unwrap();
        close(maintained.contact.normal_force_n, 0.0, 0.0, 1e-14);
        close(maintained.approach_m, dr[12] - 0.2 * c * p, 2e-12, 1e-18);
        assert_eq!(before, saved);
        assert_eq!(
            advance.contact.candidate.plastic_depth_m(),
            before.plastic_depth_m()
        );
        assert!(derivative.force_residual_n <= load_controls().maximum_force_residual_n);
    });
}

#[test]
fn finite_spring_tangent_resolves_elastic_plastic_and_detached_states() {
    with_cx(|cx, _| {
        let c = self_compliance();
        let area = 1e-6;
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let before = model.initial_state();
        let dh = vec![-2e5; 25];
        let mut dr = vec![0.0; 25];
        dr[12] = 2e-7;
        let direction = FiniteGapPlasticDriveDirection {
            position_m: 4e-7,
            stiffness_n_per_m: 2e4,
            reduced_modulus_pa: 1e8,
            hardness_pa: &dh,
            recession_m: &dr,
        };
        assert!(matches!(
            model.advance_at_load(&before, 1.0, equilibrium_controls(), cx),
            Err(FiniteGapPlasticLoadError::NonUniqueLoad { .. })
        ));
        for position_m in [3e-6, 1e-5] {
            let boundary = FiniteGapPlasticDrive {
                stiffness_n_per_m: 2e5,
                position_m,
            };
            let advance = model
                .advance_with_drive(&before, boundary, equilibrium_controls(), cx)
                .unwrap();
            let tangent = model
                .linearize_drive_advance(&before, &advance, drive_controls(), cx)
                .unwrap();
            let derivative = tangent.apply(direction, cx).unwrap();
            let plastic = position_m > 5e-6;
            let p = if plastic {
                1e6
            } else {
                boundary.stiffness_n_per_m * position_m / (area + boundary.stiffness_n_per_m * c)
            };
            let compression = area * p / boundary.stiffness_n_per_m;
            let effective = direction.position_m + compression * 0.1;
            let dp = if plastic {
                dh[12]
            } else {
                (effective - dr[12] + 0.1 * c * p) / (c + area / boundary.stiffness_n_per_m)
            };
            let da = effective - area * dp / boundary.stiffness_n_per_m;
            let du = c * (dp - 0.1 * p);
            let dz = if plastic { da - dr[12] - du } else { 0.0 };
            let increment = if plastic {
                position_m - compression - c * p
            } else {
                0.0
            };
            close(derivative.contact.pressure_pa[12], dp, 2e-12, 1e-8);
            close(derivative.approach_m, da, 2e-12, 1e-18);
            close(derivative.contact.elastic_opening_m[12], du, 2e-12, 1e-18);
            close(derivative.contact.plastic_depth_m[12], dz, 2e-12, 1e-18);
            close(
                derivative.contact.plastic_work_j,
                area * (dh[12] * increment + p * dz),
                2e-12,
                1e-18,
            );
            close(
                derivative.spring_energy_j,
                compression * (area * dp - 0.5 * compression * direction.stiffness_n_per_m),
                2e-12,
                1e-18,
            );
            assert_eq!(tangent.margins().flowing_cells, usize::from(plastic));
            assert!(derivative.force_residual_n <= drive_controls().maximum_force_residual_n);
        }
        let boundary = FiniteGapPlasticDrive {
            stiffness_n_per_m: 2e5,
            position_m: -1e-6,
        };
        let detached = model
            .advance_with_drive(&before, boundary, equilibrium_controls(), cx)
            .unwrap();
        let tangent = model
            .linearize_drive_advance(&before, &detached, drive_controls(), cx)
            .unwrap();
        let derivative = tangent.apply(direction, cx).unwrap();
        assert_eq!(tangent.margins().open_cells, 25);
        assert_eq!(derivative.approach_m, direction.position_m);
        assert_eq!(derivative.contact.pressure_pa, vec![0.0; 25]);
        assert_eq!(derivative.contact.plastic_work_j, 0.0);
        assert_eq!(derivative.spring_energy_j, 0.0);
    });
}

fn check_contact_difference(
    derivative: &FiniteGapPlasticSensitivity,
    plus: &FiniteGapPlasticAdvance,
    minus: &FiniteGapPlasticAdvance,
    epsilon: f64,
) {
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

#[test]
fn mixed_force_and_drive_tangents_match_actual_resolves_in_every_solver_lane() {
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
            let held = model.advance(&before, 3e-6, cx).unwrap();
            let force = held.receipt.normal.normal_force_n;
            let boundary = FiniteGapPlasticDrive {
                stiffness_n_per_m: 2e5,
                position_m: 3e-6 + force / 2e5,
            };
            let load = model
                .advance_at_load(&before, force, equilibrium_controls(), cx)
                .unwrap();
            let drive = model
                .advance_with_drive(&before, boundary, equilibrium_controls(), cx)
                .unwrap();
            let load_tangent = model
                .linearize_load_advance(&before, &load, load_controls(), cx)
                .unwrap();
            let drive_tangent = model
                .linearize_drive_advance(&before, &drive, drive_controls(), cx)
                .unwrap();
            assert_eq!(load_tangent.margins().elastic_cells, 1);
            assert_eq!(load_tangent.margins().flowing_cells, 1);
            assert_eq!(drive_tangent.margins().elastic_cells, 1);
            assert_eq!(drive_tangent.margins().flowing_cells, 1);
            let mut dh = vec![0.0; 25];
            dh[11] = -2e5;
            dh[13] = 3e5;
            let mut dr = vec![0.0; 25];
            dr[11] = 2e-7;
            dr[13] = -1e-7;
            let de = 2e8;
            let df = 0.1;
            let dx = 4e-7;
            let dk = 3e4;
            let load_derivative = load_tangent
                .apply(
                    FiniteGapPlasticLoadDirection {
                        target_force_n: df,
                        reduced_modulus_pa: de,
                        hardness_pa: &dh,
                        recession_m: &dr,
                    },
                    cx,
                )
                .unwrap();
            let drive_derivative = drive_tangent
                .apply(
                    FiniteGapPlasticDriveDirection {
                        position_m: dx,
                        stiffness_n_per_m: dk,
                        reduced_modulus_pa: de,
                        hardness_pa: &dh,
                        recession_m: &dr,
                    },
                    cx,
                )
                .unwrap();
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
                let loaded = owner
                    .advance_at_load(&retained, force + epsilon * df, equilibrium_controls(), cx)
                    .unwrap();
                let driven = owner
                    .advance_with_drive(
                        &retained,
                        FiniteGapPlasticDrive {
                            stiffness_n_per_m: boundary.stiffness_n_per_m + epsilon * dk,
                            position_m: boundary.position_m + epsilon * dx,
                        },
                        equilibrium_controls(),
                        cx,
                    )
                    .unwrap();
                (loaded, driven)
            };
            let epsilon = 1e-3;
            let (load_plus, drive_plus) = resolve(epsilon);
            let (load_minus, drive_minus) = resolve(-epsilon);
            check_contact_difference(
                &load_derivative.contact,
                &load_plus.contact,
                &load_minus.contact,
                epsilon,
            );
            check_contact_difference(
                &drive_derivative.contact,
                &drive_plus.contact,
                &drive_minus.contact,
                epsilon,
            );
            close(
                load_derivative.approach_m,
                (load_plus.contact.candidate.approach_m()
                    - load_minus.contact.candidate.approach_m())
                    / (2.0 * epsilon),
                2e-6,
                1e-13,
            );
            close(
                drive_derivative.approach_m,
                (drive_plus.contact.candidate.approach_m()
                    - drive_minus.contact.candidate.approach_m())
                    / (2.0 * epsilon),
                2e-6,
                1e-13,
            );
            close(
                drive_derivative.spring_energy_j,
                (drive_plus.drive.spring_energy_j - drive_minus.drive.spring_energy_j)
                    / (2.0 * epsilon),
                2e-6,
                1e-13,
            );
        }
    });
}

#[test]
fn boundary_ownership_zero_stiffness_switches_and_cancellation_are_transactional() {
    with_cx(|cx, gate| {
        let model = FiniteGapPlasticModel::new(config(), cx).unwrap();
        let before = model.initial_state();
        let saved = before.clone();
        let load = model
            .advance_at_load(&before, 0.4, equilibrium_controls(), cx)
            .unwrap();
        let boundary = FiniteGapPlasticDrive {
            stiffness_n_per_m: 2e5,
            position_m: 1e-5,
        };
        let drive = model
            .advance_with_drive(&before, boundary, equilibrium_controls(), cx)
            .unwrap();
        let mut unresolved = load_controls();
        unresolved.minimum_contact_stiffness_n_per_m = 1e20;
        assert!(matches!(
            model.linearize_load_advance(&before, &load, unresolved, cx),
            Err(FiniteGapPlasticSensitivityError::UnresolvedLoadStiffness { .. })
        ));
        let separated = model.advance(&before, -2e-6, cx).unwrap().candidate;
        let unloaded = model
            .advance_at_load(&separated, 0.0, equilibrium_controls(), cx)
            .unwrap();
        assert!(matches!(
            model.linearize_load_advance(&separated, &unloaded, load_controls(), cx),
            Err(FiniteGapPlasticSensitivityError::UnresolvedLoadStiffness {
                stiffness_n_per_m: 0.0,
                ..
            })
        ));
        let corner = model
            .advance_with_drive(
                &drive.contact.candidate,
                boundary,
                equilibrium_controls(),
                cx,
            )
            .unwrap();
        assert!(matches!(
            model.linearize_drive_advance(&drive.contact.candidate, &corner, drive_controls(), cx),
            Err(FiniteGapPlasticSensitivityError::SwitchingCell { cell: 12, .. })
        ));
        let mut wrong_load = load.clone();
        wrong_load.load.target_force_n = 10.0;
        assert!(
            model
                .linearize_load_advance(&before, &wrong_load, load_controls(), cx)
                .is_err()
        );
        wrong_load = load.clone();
        wrong_load.load.predecessor_state_identity = separated.identity();
        assert!(
            model
                .linearize_load_advance(&before, &wrong_load, load_controls(), cx)
                .is_err()
        );
        let mut wrong_drive = drive.clone();
        wrong_drive.drive.boundary.stiffness_n_per_m = 0.0;
        assert!(
            model
                .linearize_drive_advance(&before, &wrong_drive, drive_controls(), cx)
                .is_err()
        );
        wrong_drive = drive.clone();
        wrong_drive.contact.receipt.predecessor_identity = separated.identity();
        assert!(
            model
                .linearize_drive_advance(&before, &wrong_drive, drive_controls(), cx)
                .is_err()
        );
        let tangent = model
            .linearize_drive_advance(&before, &drive, drive_controls(), cx)
            .unwrap();
        let zero = vec![0.0; 25];
        let direction = FiniteGapPlasticDriveDirection {
            position_m: 0.0,
            stiffness_n_per_m: 0.0,
            reduced_modulus_pa: 0.0,
            hardness_pa: &zero,
            recession_m: &zero,
        };
        assert!(
            tangent
                .apply(
                    FiniteGapPlasticDriveDirection {
                        stiffness_n_per_m: f64::NAN,
                        ..direction
                    },
                    cx
                )
                .is_err()
        );
        assert!(
            tangent
                .apply(
                    FiniteGapPlasticDriveDirection {
                        hardness_pa: &zero[..24],
                        ..direction
                    },
                    cx
                )
                .is_err()
        );
        gate.request();
        assert!(tangent.apply(direction, cx).is_err());
        assert_eq!(before, saved);
    });
}
