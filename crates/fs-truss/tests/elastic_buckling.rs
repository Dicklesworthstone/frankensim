//! G1/G3 global truss buckling with fixed thermal/settlement prestress.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_truss::elastic::buckling::{BucklingControls, BucklingError, BucklingLoadFamily, solve};
use fs_truss::elastic::{ElasticControls, ElasticError, ElasticImposedState};
use fs_truss::{GroundLimits, GroundStructure, LayoutCase, TrussConstructionError};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 1,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-8 * expected.abs().max(1e-15),
        "{actual:e} != {expected:e}",
    );
}

fn braced_node(cx: &Cx<'_>, angle: f64, reference_force: f64) -> (GroundStructure, LayoutCase) {
    let (s, c) = angle.sin_cos();
    let rotate = |p: [f64; 2]| [c * p[0] - s * p[1], s * p[0] + c * p[1]];
    let ground = GroundStructure::try_from_parts(
        &[[-2.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, -1.0]].map(rotate),
        &[(0, 1), (1, 2), (1, 3)],
        &[2.0, 1.0, 1.0],
        GroundLimits::default(),
        cx,
    )
    .unwrap();
    let case = LayoutCase::try_new(
        vec![[true; 2], [false; 2], [true; 2], [true; 2]],
        vec![[0.0; 2], rotate([reference_force, 0.0]), [0.0; 2], [0.0; 2]],
        4,
    )
    .unwrap();
    (ground, case)
}

#[test]
fn g1_fixed_thermal_settlement_and_dead_loads_shift_the_global_critical_factor() {
    with_cx(|cx, _| {
        // Both horizontal EA/L = k=1e6; vertical brace EA/L=100.
        // At fixed baseline N0=k/2*(settlement-3*epsilon), a fixed
        // horizontal force F adds [+F/2,-F/2] to the horizontal members.
        // K_base,yy=100+1.5*N0-F/4, K_reference,yy=-P/4.
        for (epsilon, settlement, fixed_force) in [
            (0.0, 0.0, 0.0),
            (1e-5, 0.0, 0.0),
            (1e-5, 2e-5, 20.0),
            (1e-5, -2e-5, 0.0),
            (-1e-5, 0.0, 0.0),
        ] {
            for reference_force in [1.0, 7.0] {
                let (ground, case) = braced_node(cx, 0.0, reference_force);
                let strains = [epsilon, epsilon, 0.0];
                let motion = [[0.0; 2], [0.0; 2], [settlement, 0.0], [0.0; 2]];
                let fixed_loads = [[0.0; 2], [fixed_force, 0.0], [0.0; 2], [0.0; 2]];
                let result = solve(
                    &ground,
                    &case,
                    &[1.0; 3],
                    &[2e6, 1e6, 100.0],
                    BucklingLoadFamily {
                        imposed: ElasticImposedState {
                            member_stress_free_strain: Some(&strains),
                            prescribed_displacements_m: Some(&motion),
                        },
                        fixed_loads_n: Some(&fixed_loads),
                    },
                    BucklingControls::default(),
                    cx,
                )
                .unwrap();
                let n0 = 0.5e6 * (settlement - 3.0 * epsilon);
                let expected = 4.0 * (100.0 + 1.5 * n0 - fixed_force / 4.0) / reference_force;
                assert_eq!(result.buckling.loads.len(), 1);
                near(result.buckling.loads[0], expected);
                near(
                    result.reference_load.member_force_n[0],
                    reference_force / 2.0,
                );
                near(
                    result.reference_load.member_force_n[1],
                    -reference_force / 2.0,
                );
                assert!((result.baseline.member_force_n[0] - n0 - fixed_force / 2.0).abs() < 1e-10);
                assert!((result.baseline.member_force_n[1] - n0 + fixed_force / 2.0).abs() < 1e-10);
                assert_eq!(result.baseline.displacements_m[2], motion[2]);
                assert_eq!(result.reference_load.displacements_m[2], [0.0; 2]);
                assert_eq!(result.buckling.dof_map[2], Some(0));
                assert_eq!(result.buckling.dof_map[3], Some(1));
                assert_eq!(result.buckling.modes[0].len(), 2);
                assert!(result.buckling.relative_residuals[0] <= 1e-8);
            }
        }
    });
}

#[test]
fn g3_rotation_preserves_the_biased_load_factor_and_rotates_the_global_mode() {
    with_cx(|cx, _| {
        for angle in [0.0_f64, 0.37, 1.1] {
            let (ground, case) = braced_node(cx, angle, 1.0);
            let (s, c) = angle.sin_cos();
            let motion = [[0.0; 2], [0.0; 2], [2e-5 * c, 2e-5 * s], [0.0; 2]];
            let fixed_loads = [[0.0; 2], [20.0 * c, 20.0 * s], [0.0; 2], [0.0; 2]];
            let result = solve(
                &ground,
                &case,
                &[1.0; 3],
                &[2e6, 1e6, 100.0],
                BucklingLoadFamily {
                    imposed: ElasticImposedState {
                        member_stress_free_strain: Some(&[1e-5, 1e-5, 0.0]),
                        prescribed_displacements_m: Some(&motion),
                    },
                    fixed_loads_n: Some(&fixed_loads),
                },
                BucklingControls::default(),
                cx,
            )
            .unwrap();
            near(result.buckling.loads[0], 350.0);
            let mode = &result.buckling.modes[0];
            let amplitude = mode[0].hypot(mode[1]);
            assert!(amplitude > 0.0);
            assert!((c * mode[0] + s * mode[1]).abs() < 1e-8 * amplitude);
            assert!(result.buckling.relative_residuals[0] <= 1e-8);
        }
    });
}

#[test]
fn g1_tensile_and_zero_reference_loads_return_no_positive_buckling_root() {
    with_cx(|cx, _| {
        let ground = GroundStructure::try_from_parts(
            &[[-1.0, 0.0], [0.0, -1.0], [0.0, 0.0]],
            &[(0, 2), (1, 2)],
            &[1.0, 1.0],
            GroundLimits::default(),
            cx,
        )
        .unwrap();
        for force in [0.0, 1.0] {
            let case = LayoutCase::try_new(
                vec![[true; 2], [true; 2], [false; 2]],
                vec![[0.0; 2], [0.0; 2], [force; 2]],
                3,
            )
            .unwrap();
            let result = solve(
                &ground,
                &case,
                &[1.0; 2],
                &[1e6; 2],
                BucklingLoadFamily::default(),
                BucklingControls::default(),
                cx,
            )
            .unwrap();
            assert!(result.buckling.loads.is_empty());
            assert!(result.buckling.modes.is_empty());
            for member_force in result.reference_load.member_force_n {
                assert!((member_force - force).abs() < 1e-12);
            }
        }
    });
}

#[test]
#[allow(clippy::too_many_lines)] // Independent refusals on one physical fixture.
fn g4_mechanisms_unstable_baseline_invalid_data_work_and_cancellation_refuse() {
    with_cx(|cx, gate| {
        let (ground, case) = braced_node(cx, 0.0, 1.0);
        let young = [2e6, 1e6, 100.0];
        let controls = BucklingControls::default();
        let family = BucklingLoadFamily::default();
        // The ordinary axial solve can admit this unloaded transverse mechanism.
        assert!(
            fs_truss::elastic::solve(
                &ground,
                &case,
                &[1.0, 1.0, 0.0],
                &young,
                ElasticControls::default(),
                cx,
            )
            .is_ok()
        );
        // Stability retains that free coordinate and refuses its singular base.
        assert!(matches!(
            solve(
                &ground,
                &case,
                &[1.0, 1.0, 0.0],
                &young,
                family,
                controls,
                cx
            ),
            Err(BucklingError::Stability(
                fs_solid::SolidError::SolveFailed { .. }
            ))
        ));
        assert!(matches!(
            solve(
                &ground,
                &case,
                &[1.0; 3],
                &young,
                BucklingLoadFamily {
                    imposed: ElasticImposedState {
                        member_stress_free_strain: Some(&[1e-4, 1e-4, 0.0]),
                        ..ElasticImposedState::default()
                    },
                    ..family
                },
                controls,
                cx,
            ),
            Err(BucklingError::Stability(
                fs_solid::SolidError::SolveFailed { .. }
            ))
        ));
        assert!(
            solve(
                &ground,
                &case,
                &[1.0; 3],
                &young,
                BucklingLoadFamily {
                    fixed_loads_n: Some(&[[0.0; 2]; 3]),
                    ..family
                },
                controls,
                cx,
            )
            .is_err()
        );
        assert!(
            solve(
                &ground,
                &case,
                &[1.0; 3],
                &young,
                family,
                BucklingControls {
                    modes: 0,
                    ..controls
                },
                cx,
            )
            .is_err()
        );
        assert!(matches!(
            solve(
                &ground,
                &case,
                &[1.0; 3],
                &young,
                family,
                BucklingControls {
                    elastic: ElasticControls {
                        maximum_iterations: 0,
                        ..controls.elastic
                    },
                    ..controls
                },
                cx,
            ),
            Err(BucklingError::Elastic(ElasticError::NotConverged { .. }))
        ));
        assert!(matches!(
            solve(
                &ground,
                &case,
                &[1.0; 3],
                &young,
                family,
                BucklingControls {
                    maximum_staged_entries: 1,
                    ..controls
                },
                cx,
            ),
            Err(BucklingError::Elastic(ElasticError::Construction(
                TrussConstructionError::WorkBudget { .. }
            )))
        ));
        gate.request();
        assert!(matches!(
            solve(&ground, &case, &[1.0; 3], &young, family, controls, cx),
            Err(BucklingError::Elastic(ElasticError::Construction(
                TrussConstructionError::Cancelled { .. }
            )))
        ));
    });
}
