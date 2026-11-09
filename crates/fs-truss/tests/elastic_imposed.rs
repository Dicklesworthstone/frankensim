//! Analytic thermal/misfit strain and settlement through the native solve.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_truss::elastic::{
    ElasticControls, ElasticError, ElasticImposedState, solve, solve_with_imposed,
};
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
        (actual - expected).abs() <= 2e-11 * actual.abs().max(expected.abs()).max(1e-15),
        "{actual:e} != {expected:e}",
    );
}

fn bar(cx: &Cx<'_>) -> GroundStructure {
    GroundStructure::try_from_parts(
        &[[0.0, 0.0], [2.0, 0.0]],
        &[(0, 1)],
        &[2.0],
        GroundLimits::default(),
        cx,
    )
    .unwrap()
}

fn series(cx: &Cx<'_>) -> GroundStructure {
    GroundStructure::try_from_parts(
        &[[0.0, 0.0], [2.0, 0.0], [3.0, 0.0]],
        &[(0, 1), (1, 2)],
        &[2.0, 1.0],
        GroundLimits::default(),
        cx,
    )
    .unwrap()
}

#[test]
fn free_thermal_expansion_and_rigid_translation_store_no_elastic_energy() {
    with_cx(|cx, _| {
        let ground = bar(cx);
        let case =
            LayoutCase::try_new(vec![[true, true], [false, true]], vec![[0.0; 2]; 2], 2).unwrap();
        for translation in [0.0, 3e-5] {
            let motion = [[translation, 0.0], [0.0; 2]];
            let state = ElasticImposedState {
                member_stress_free_strain: Some(&[1e-5]),
                prescribed_displacements_m: Some(&motion),
            };
            let response = solve_with_imposed(
                &ground,
                &case,
                &[2e-4],
                &[200e9],
                state,
                ElasticControls::default(),
                cx,
            )
            .unwrap();
            assert_eq!(response.displacements_m[0], motion[0]);
            near(response.displacements_m[1][0], translation + 2e-5);
            near(response.member_strain[0], 1e-5);
            assert!(response.member_elastic_strain[0].abs() < 1e-18);
            assert!(response.member_force_n[0].abs() < 1e-9);
            assert!(response.support_reactions_n[0][0].abs() < 1e-9);
            assert!(response.strain_energy_j < 1e-20);
            assert_eq!(response.compliance_j, 0.0);
            assert!(response.relative_equilibrium_residual < 1e-12);
        }
    });
}

#[test]
fn rotated_bar_prescribed_rigid_translation_has_exactly_zero_stress_and_energy() {
    with_cx(|cx, _| {
        let ground = GroundStructure::try_from_parts(
            &[[0.0, 0.0], [3.0, 5.0]],
            &[(0, 1)],
            &[3.0_f64.hypot(5.0)],
            GroundLimits::default(),
            cx,
        )
        .unwrap();
        let case = LayoutCase::try_new(vec![[true; 2]; 2], vec![[0.0; 2]; 2], 2).unwrap();
        let motion = [[1e6; 2]; 2];
        let response = solve_with_imposed(
            &ground,
            &case,
            &[2e-4],
            &[200e9],
            ElasticImposedState {
                prescribed_displacements_m: Some(&motion),
                ..ElasticImposedState::default()
            },
            ElasticControls::default(),
            cx,
        )
        .unwrap();
        // Large common translation is a numerical invariance check: there is
        // no relative motion or deformation, despite non-axis-aligned members.
        assert_eq!(response.displacements_m, motion.to_vec());
        assert_eq!(response.member_strain, vec![0.0]);
        assert_eq!(response.member_elastic_strain, vec![0.0]);
        assert_eq!(response.member_force_n, vec![0.0]);
        assert_eq!(response.member_stress_pa, vec![0.0]);
        assert_eq!(response.support_reactions_n, vec![[0.0; 2]; 2]);
        assert_eq!(response.strain_energy_j, 0.0);
        assert_eq!(response.compliance_j, 0.0);
        assert_eq!(response.iterations, 0);
    });
}

#[test]
fn fully_restrained_thermal_expansion_recovers_stress_reactions_and_energy() {
    with_cx(|cx, _| {
        let ground = bar(cx);
        let case = LayoutCase::try_new(vec![[true; 2]; 2], vec![[0.0; 2]; 2], 2).unwrap();
        for sign in [-1.0, 1.0] {
            let strains = [sign * 1e-5];
            let response = solve_with_imposed(
                &ground,
                &case,
                &[2e-4],
                &[200e9],
                ElasticImposedState {
                    member_stress_free_strain: Some(&strains),
                    ..ElasticImposedState::default()
                },
                ElasticControls::default(),
                cx,
            )
            .unwrap();
            assert_eq!(response.displacements_m, vec![[0.0; 2]; 2]);
            assert_eq!(response.member_strain[0], 0.0);
            near(response.member_elastic_strain[0], -sign * 1e-5);
            near(response.member_force_n[0], -sign * 400.0);
            near(response.member_stress_pa[0], -sign * 2e6);
            near(response.support_reactions_n[0][0], sign * 400.0);
            near(response.support_reactions_n[1][0], -sign * 400.0);
            near(response.strain_energy_j, 0.004);
            assert_eq!(response.compliance_j, 0.0);
            assert_eq!(response.iterations, 0);
        }
    });
}

#[test]
fn heterogeneous_series_keeps_one_force_and_distinct_total_and_elastic_strains() {
    with_cx(|cx, _| {
        let ground = series(cx);
        let case = LayoutCase::try_new(
            vec![[true, true], [false, true], [false, true]],
            vec![[0.0; 2], [0.0; 2], [12.0, 0.0]],
            3,
        )
        .unwrap();
        let response = solve_with_imposed(
            &ground,
            &case,
            &[2.0, 3.0],
            &[4e6, 2e6],
            ElasticImposedState {
                member_stress_free_strain: Some(&[7e-6, -3e-6]),
                ..ElasticImposedState::default()
            },
            ElasticControls::default(),
            cx,
        )
        .unwrap();
        // Extensions are L*epsilon_0 + P*L/(EA): 17e-6 and -1e-6.
        near(response.displacements_m[1][0], 17e-6);
        near(response.displacements_m[2][0], 16e-6);
        near(response.member_strain[0], 8.5e-6);
        near(response.member_strain[1], -1e-6);
        near(response.member_elastic_strain[0], 1.5e-6);
        near(response.member_elastic_strain[1], 2e-6);
        near(response.member_force_n[0], 12.0);
        near(response.member_force_n[1], 12.0);
        near(response.support_reactions_n[0][0], -12.0);
        near(response.strain_energy_j, 30e-6);
        near(response.compliance_j, 192e-6);
        assert!(response.relative_equilibrium_residual < 1e-12);
    });
}

#[test]
fn settlement_external_loads_and_misfit_satisfy_equilibrium_and_full_work_identity() {
    with_cx(|cx, _| {
        let ground = series(cx);
        let case = LayoutCase::try_new(
            vec![[true, true], [false, true], [true, true]],
            vec![[2.0, 0.0], [10.0, 0.0], [-3.0, 0.0]],
            3,
        )
        .unwrap();
        let motion = [[-2e-6, 0.0], [0.0; 2], [8e-6, 0.0]];
        let strains = [1e-6, -2e-6];
        let response = solve_with_imposed(
            &ground,
            &case,
            &[2.0, 3.0],
            &[4e6, 2e6],
            ElasticImposedState {
                member_stress_free_strain: Some(&strains),
                prescribed_displacements_m: Some(&motion),
            },
            ElasticControls::default(),
            cx,
        )
        .unwrap();
        // N1 = 4e6*u1, N2 = 6e6*(10e-6-u1), N1-N2 = 10.
        assert_eq!(response.displacements_m[0], motion[0]);
        assert_eq!(response.displacements_m[2], motion[2]);
        near(response.displacements_m[1][0], 7e-6);
        near(response.member_force_n[0], 28.0);
        near(response.member_force_n[1], 18.0);
        near(response.member_elastic_strain[0], 3.5e-6);
        near(response.member_elastic_strain[1], 3e-6);
        near(response.support_reactions_n[0][0], -30.0);
        near(response.support_reactions_n[2][0], 21.0);
        near(response.strain_energy_j, 125e-6);
        near(response.compliance_j, 42e-6);
        let support_work: f64 = response
            .support_reactions_n
            .iter()
            .zip(motion)
            .map(|(reaction, displacement)| reaction[0] * displacement[0])
            .sum();
        let eigenstrain_work: f64 = response
            .member_force_n
            .iter()
            .zip(ground.lengths())
            .zip(strains)
            .map(|((&force, &length), strain)| force * length * strain)
            .sum();
        near(
            2.0 * response.strain_energy_j,
            response.compliance_j + support_work - eigenstrain_work,
        );
        let resultant: f64 = response
            .support_reactions_n
            .iter()
            .zip(case.loads())
            .map(|(reaction, load)| reaction[0] + load[0])
            .sum();
        assert!(resultant.abs() < 1e-11);
        assert!(response.relative_equilibrium_residual < 1e-12);
    });
}

#[test]
fn balancing_external_and_settlement_loads_do_not_break_residual_normalization() {
    with_cx(|cx, _| {
        let ground = bar(cx);
        let case = LayoutCase::try_new(
            vec![[true, true], [false, true]],
            vec![[0.0; 2], [-400.0, 0.0]],
            2,
        )
        .unwrap();
        let response = solve_with_imposed(
            &ground,
            &case,
            &[2e-4],
            &[200e9],
            ElasticImposedState {
                prescribed_displacements_m: Some(&[[2e-5, 0.0], [0.0; 2]]),
                ..ElasticImposedState::default()
            },
            ElasticControls::default(),
            cx,
        )
        .unwrap();
        assert!(response.displacements_m[1][0].abs() < 1e-18);
        near(response.member_force_n[0], -400.0);
        near(response.support_reactions_n[0][0], 400.0);
        near(response.strain_energy_j, 0.004);
        assert!(response.relative_equilibrium_residual < 1e-12);
    });
}

#[test]
fn zero_imposed_state_preserves_original_results_and_removed_members_store_no_stress() {
    with_cx(|cx, _| {
        let ground = bar(cx);
        let case = LayoutCase::try_new(
            vec![[true, true], [false, true]],
            vec![[0.0; 2], [12.0, 0.0]],
            2,
        )
        .unwrap();
        let original = solve(
            &ground,
            &case,
            &[2.0],
            &[4e6],
            ElasticControls::default(),
            cx,
        )
        .unwrap();
        let explicit = solve_with_imposed(
            &ground,
            &case,
            &[2.0],
            &[4e6],
            ElasticImposedState {
                member_stress_free_strain: Some(&[0.0]),
                prescribed_displacements_m: Some(&[[0.0; 2]; 2]),
            },
            ElasticControls::default(),
            cx,
        )
        .unwrap();
        assert_eq!(original.displacements_m, explicit.displacements_m);
        assert_eq!(original.member_strain, explicit.member_strain);
        assert_eq!(original.member_strain, original.member_elastic_strain);
        assert_eq!(original.member_force_n, explicit.member_force_n);
        assert_eq!(original.member_stress_pa, explicit.member_stress_pa);
        assert_eq!(original.support_reactions_n, explicit.support_reactions_n);
        assert_eq!(original.strain_energy_j, explicit.strain_energy_j);
        assert_eq!(original.compliance_j, explicit.compliance_j);
        assert_eq!(original.iterations, explicit.iterations);

        let unloaded =
            LayoutCase::try_new(vec![[true, true], [false, true]], vec![[0.0; 2]; 2], 2).unwrap();
        let removed = solve_with_imposed(
            &ground,
            &unloaded,
            &[0.0],
            &[4e6],
            ElasticImposedState {
                member_stress_free_strain: Some(&[f64::MAX]),
                prescribed_displacements_m: Some(&[[3e-5, 0.0], [0.0; 2]]),
            },
            ElasticControls::default(),
            cx,
        )
        .unwrap();
        assert_eq!(removed.displacements_m, vec![[3e-5, 0.0], [0.0; 2]]);
        assert_eq!(removed.member_strain, vec![0.0]);
        assert_eq!(removed.member_elastic_strain, vec![0.0]);
        assert_eq!(removed.member_force_n, vec![0.0]);
        assert_eq!(removed.strain_energy_j, 0.0);
    });
}

#[test]
fn malformed_imposed_state_nonfinite_lifting_and_cancellation_refuse() {
    with_cx(|cx, gate| {
        let ground = bar(cx);
        let case =
            LayoutCase::try_new(vec![[true, true], [false, true]], vec![[0.0; 2]; 2], 2).unwrap();
        let run = |state| {
            solve_with_imposed(
                &ground,
                &case,
                &[2e-4],
                &[200e9],
                state,
                ElasticControls::default(),
                cx,
            )
        };
        for strains in [&[0.0, 0.0][..], &[f64::NAN][..], &[f64::MAX][..]] {
            assert!(
                run(ElasticImposedState {
                    member_stress_free_strain: Some(strains),
                    ..ElasticImposedState::default()
                })
                .is_err()
            );
        }
        for motion in [
            &[[0.0; 2]][..],
            &[[0.0; 2], [1e-6, 0.0]][..],
            &[[f64::NAN, 0.0], [0.0; 2]][..],
            &[[f64::MAX, 0.0], [0.0; 2]][..],
        ] {
            assert!(
                run(ElasticImposedState {
                    prescribed_displacements_m: Some(motion),
                    ..ElasticImposedState::default()
                })
                .is_err()
            );
        }
        gate.request();
        assert!(matches!(
            run(ElasticImposedState::default()),
            Err(ElasticError::Construction(
                TrussConstructionError::Cancelled { .. }
            ))
        ));
    });
}
