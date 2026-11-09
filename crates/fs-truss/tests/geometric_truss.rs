//! G1/G3 objective network equilibrium, exact tangent and nonlinear folds.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::SolidError;
use fs_solid::continuation::{ArcSettings, PathEvent};
use fs_solver::nonlinear::NewtonKrylovConfig;
use fs_truss::elastic::ElasticImposedState;
use fs_truss::geometric::{GeometricLimits, GeometricLoading, GeometricScales, GeometricTruss};
use fs_truss::{GroundLimits, GroundStructure, LayoutCase};

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
fn limits() -> GeometricLimits {
    GeometricLimits {
        maximum_absolute_elastic_strain: 0.2,
        maximum_free_dofs: 100,
        maximum_staged_entries: 10000,
        maximum_payload_bytes: 16000000,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}"
    );
}
fn ground(nodes: &[[f64; 2]], members: &[(usize, usize)], cx: &Cx<'_>) -> GroundStructure {
    let lengths: Vec<_> = members
        .iter()
        .map(|&(a, b)| (nodes[b][0] - nodes[a][0]).hypot(nodes[b][1] - nodes[a][1]))
        .collect();
    GroundStructure::try_from_parts(nodes, members, &lengths, GroundLimits::default(), cx).unwrap()
}

#[test]
fn g1_native_equilibrium_keeps_thermal_strain_support_ramp_and_supported_loads() {
    with_cx(|cx, _| {
        let g = ground(&[[0.0, 0.0], [2.0, 0.0]], &[(0, 1)], cx);
        let case = LayoutCase::try_new(
            vec![[true; 2], [false, true]],
            vec![[7.0, -1.0], [10.0, 5.0]],
            2,
        )
        .unwrap();
        let fixed_motion = [[0.02, 0.01], [0.0, 0.01]];
        let ramp = [[0.003, 0.002], [0.0, 0.002]];
        let fixed_loads = [[3.0, 4.0], [5.0, 6.0]];
        let truss = GeometricTruss::try_new(
            &g,
            &case,
            &[2.0],
            &[1000.0],
            GeometricLoading {
                fixed: ElasticImposedState {
                    member_stress_free_strain: Some(&[0.01]),
                    prescribed_displacements_m: Some(&fixed_motion),
                },
                fixed_loads_n: Some(&fixed_loads),
                support_motion_per_factor_m: Some(&ramp),
            },
            GeometricScales {
                displacement_m: 0.1,
                force_n: 100.0,
            },
            limits(),
            cx,
        )
        .unwrap();
        let support_x = fixed_motion[0][0] + 2.0 * ramp[0][0];
        let support_y = fixed_motion[0][1] + 2.0 * ramp[0][1];
        let seed = [[support_x, support_y], [0.0, support_y]];
        let response = truss
            .equilibrate(&seed, 2.0, NewtonKrylovConfig::default(), 12, cx)
            .unwrap();
        near(response.displacements_m[1][0], 0.071, 2e-11);
        near(response.member_engineering_strain[0], 0.0225, 2e-11);
        near(response.member_elastic_strain[0], 0.0125, 2e-11);
        near(response.member_force_n[0], 25.0, 2e-11);
        near(response.member_nominal_stress_pa[0], 12.5, 2e-11);
        near(response.strain_energy_j, 0.3125, 2e-11);
        near(response.support_reactions_n[0][0], -42.0, 2e-11);
        near(response.support_reactions_n[0][1], -2.0, 2e-11);
        near(response.support_reactions_n[1][1], -16.0, 2e-11);
        assert!(response.free_force_residual_n[0].abs() < 1e-9);
        assert!(
            truss
                .equilibrate(&seed, 2.0, NewtonKrylovConfig::default(), 0, cx)
                .is_err()
        );
        let path = truss
            .start_from(&response.displacements_m, 2.0, &ArcSettings::default(), cx)
            .unwrap();
        assert_eq!(path.lambda, 2.0);
        assert_eq!(truss.free_dofs(), &[2]);
    });
}

#[test]
fn g3_current_chord_tangent_load_column_and_energy_match_independent_differences() {
    with_cx(|cx, _| {
        let g = ground(
            &[[-1.0, 0.0], [1.0, 0.0], [0.2, 0.6]],
            &[(0, 2), (1, 2)],
            cx,
        );
        let case = LayoutCase::try_new(
            vec![[true; 2], [true; 2], [false; 2]],
            vec![[2.0, 3.0], [0.0; 2], [-7.0, -11.0]],
            3,
        )
        .unwrap();
        let fixed_motion = [[0.004, -0.002], [-0.003, 0.001], [0.0; 2]];
        let ramp = [[0.002, 0.003], [-0.001, 0.004], [0.0; 2]];
        let fixed_loads = [[5.0, -7.0], [0.0; 2], [3.0, -4.0]];
        let truss = GeometricTruss::try_new(
            &g,
            &case,
            &[1.2, 0.8],
            &[1200.0, 2100.0],
            GeometricLoading {
                fixed: ElasticImposedState {
                    member_stress_free_strain: Some(&[0.01, -0.005]),
                    prescribed_displacements_m: Some(&fixed_motion),
                },
                fixed_loads_n: Some(&fixed_loads),
                support_motion_per_factor_m: Some(&ramp),
            },
            GeometricScales {
                displacement_m: 0.2,
                force_n: 10.0,
            },
            limits(),
            cx,
        )
        .unwrap();
        let q = [0.03, -0.07];
        let lambda = 0.8;
        let h = 2e-5;
        let k = truss.tangent(&q, lambda, cx).unwrap();
        let f = truss.load_vector_at(&q, lambda, cx).unwrap();
        let response = truss.response(&q, lambda, cx).unwrap();
        for j in 0..2 {
            let sample = |s: f64| {
                let mut x = q;
                x[j] += s * h;
                truss.response(&x, lambda, cx).unwrap()
            };
            let [mm, m, p, pp] = [-2.0, -1.0, 1.0, 2.0].map(sample);
            for i in 0..2 {
                let derivative = (mm.free_force_residual_n[i] - 8.0 * m.free_force_residual_n[i]
                    + 8.0 * p.free_force_residual_n[i]
                    - pp.free_force_residual_n[i])
                    / (12.0 * h * 10.0);
                near(k.get(i, j), derivative, 3e-8);
                near(k.get(i, j), k.get(j, i), 1e-14);
            }
            let energy_derivative = (mm.strain_energy_j - 8.0 * m.strain_energy_j
                + 8.0 * p.strain_energy_j
                - pp.strain_energy_j)
                / (12.0 * h);
            let internal =
                response.free_force_residual_n[j] + fixed_loads[2][j] + lambda * case.loads()[2][j];
            near(energy_derivative, 0.2 * internal, 3e-8);
        }
        let r: Vec<_> = [-2.0, -1.0, 1.0, 2.0]
            .iter()
            .map(|s| truss.residual(&q, lambda + s * h, cx).unwrap())
            .collect();
        for i in 0..2 {
            near(
                -f[i],
                (r[0][i] - 8.0 * r[1][i] + 8.0 * r[2][i] - r[3][i]) / (12.0 * h),
                3e-8,
            );
            let sum = (0..3)
                .map(|node| {
                    response.support_reactions_n[node][i]
                        + fixed_loads[node][i]
                        + lambda * case.loads()[node][i]
                })
                .sum::<f64>()
                + response.free_force_residual_n[i];
            near(sum, 0.0, 1e-10);
        }
        assert!(k.get(0, 1).abs() > 1.0); // This probes a genuinely coupled tangent.
    });
}

#[test]
fn g3_finite_rigid_rotations_and_free_thermal_expansion_remain_stress_free() {
    with_cx(|cx, _| {
        let nodes = [[-1.0, 0.0], [1.0, 0.0], [0.2, 0.6]];
        let g = ground(&nodes, &[(0, 2), (1, 2)], cx);
        let case =
            LayoutCase::try_new(vec![[true; 2], [true; 2], [false; 2]], vec![[0.0; 2]; 3], 3)
                .unwrap();
        for angle in [0.0_f64, 0.4, 1.7, 3.0] {
            for expansion in [0.0, 0.015] {
                let (s, c) = angle.sin_cos();
                let u: Vec<_> = nodes
                    .iter()
                    .map(|p| {
                        [
                            (1.0 + expansion) * (c * p[0] - s * p[1]) - p[0] + 0.3,
                            (1.0 + expansion) * (s * p[0] + c * p[1]) - p[1] - 0.2,
                        ]
                    })
                    .collect();
                let supports = [u[0], u[1], [0.0; 2]];
                let eigenstrain = [expansion; 2];
                let truss = GeometricTruss::try_new(
                    &g,
                    &case,
                    &[1.0; 2],
                    &[2000.0; 2],
                    GeometricLoading {
                        fixed: ElasticImposedState {
                            member_stress_free_strain: Some(&eigenstrain),
                            prescribed_displacements_m: Some(&supports),
                        },
                        ..GeometricLoading::default()
                    },
                    GeometricScales {
                        displacement_m: 1.0,
                        force_n: 1.0,
                    },
                    limits(),
                    cx,
                )
                .unwrap();
                let q = truss.free_coordinates(&u, 0.0, cx).unwrap();
                let response = truss.response(&q, 0.0, cx).unwrap();
                assert!(response.member_force_n.iter().all(|n| n.abs() < 2e-12));
                assert!(response.strain_energy_j < 1e-25);
            }
        }
    });
}

#[test]
fn g1_native_arclength_tracks_both_shallow_truss_folds_and_resumes_exactly() {
    with_cx(|cx, _| {
        let height = 0.25_f64;
        let length = 1.0_f64.hypot(height);
        let g = ground(
            &[[-1.0, 0.0], [1.0, 0.0], [0.0, height]],
            &[(0, 2), (1, 2)],
            cx,
        );
        let case = LayoutCase::try_new(
            vec![[true; 2], [true; 2], [true, false]],
            vec![[0.0; 2], [0.0; 2], [0.0, -10.0]],
            3,
        )
        .unwrap();
        let truss = GeometricTruss::try_new(
            &g,
            &case,
            &[1.0; 2],
            &[1000.0; 2],
            GeometricLoading::default(),
            GeometricScales {
                displacement_m: height,
                force_n: 10.0,
            },
            limits(),
            cx,
        )
        .unwrap();
        let settings = ArcSettings {
            ds: 0.03,
            ds_min: 1e-6,
            ds_max: 0.06,
            max_corrector_iters: 12,
            tol: 1e-10,
        };
        let start = truss
            .start_from(&[[0.0; 2]; 3], 0.0, &settings, cx)
            .unwrap();
        let one = truss.advance(&start, &settings, 12, cx).unwrap();
        let split = truss
            .advance(
                &truss.advance(&start, &settings, 5, cx).unwrap(),
                &settings,
                7,
                cx,
            )
            .unwrap();
        assert_eq!(one.u, split.u);
        assert_eq!(one.lambda.to_bits(), split.lambda.to_bits());
        assert_eq!(one.tangent, split.tangent);
        assert_eq!(one.trace, split.trace);
        assert_eq!(one.events, split.events);
        let mut state = start;
        let mut maximum = 0.0_f64;
        let mut minimum = 0.0_f64;
        for _ in 0..160 {
            state = truss.advance(&state, &settings, 1, cx).unwrap();
            let response = truss.response(&state.u, state.lambda, cx).unwrap();
            let y = height + response.displacements_m[2][1];
            let current_length = 1.0_f64.hypot(y);
            let analytic =
                -2.0 * 1000.0 * (current_length - length) / length * y / current_length / 10.0;
            near(state.lambda, analytic, 2e-9);
            assert!(response.free_force_residual_n[0].abs() < 1.1e-9);
            maximum = maximum.max(state.lambda);
            minimum = minimum.min(state.lambda);
            if state.u[0] < -2.05 {
                break;
            }
        }
        assert!(state.u[0] < -2.05);
        assert_eq!(
            state
                .events
                .iter()
                .filter(|e| matches!(e, PathEvent::LimitPoint { .. }))
                .count(),
            2
        );
        let critical_length = length.cbrt();
        let critical_y = (critical_length * critical_length - 1.0).sqrt();
        let critical = 2.0 * 1000.0 * (length - critical_length) / length * critical_y
            / critical_length
            / 10.0;
        near(maximum, critical, 2e-3);
        near(minimum, -critical, 2e-3);
    });
}

#[test]
fn g4_domains_resources_unbalanced_seed_and_cancellation_refuse() {
    with_cx(|cx, gate| {
        let g = ground(&[[0.0, 0.0], [1.0, 0.0]], &[(0, 1)], cx);
        let case = LayoutCase::try_new(
            vec![[true; 2], [false, true]],
            vec![[0.0; 2], [1.0, 0.0]],
            2,
        )
        .unwrap();
        let scales = GeometricScales {
            displacement_m: 1.0,
            force_n: 1.0,
        };
        let truss = GeometricTruss::try_new(
            &g,
            &case,
            &[1.0],
            &[1000.0],
            GeometricLoading::default(),
            scales,
            limits(),
            cx,
        )
        .unwrap();
        assert!(truss.response(&[-1.0], 0.0, cx).is_err());
        assert!(truss.response(&[0.21], 0.0, cx).is_err());
        assert!(truss.response(&[f64::NAN], 0.0, cx).is_err());
        assert!(
            truss
                .start_from(&[[0.0; 2]; 2], 1.0, &ArcSettings::default(), cx)
                .is_err()
        );
        assert!(
            GeometricTruss::try_new(
                &g,
                &case,
                &[1.0],
                &[1000.0],
                GeometricLoading::default(),
                scales,
                GeometricLimits {
                    maximum_payload_bytes: 1,
                    ..limits()
                },
                cx
            )
            .is_err()
        );
        assert!(
            GeometricTruss::try_new(
                &g,
                &case,
                &[1.0],
                &[1000.0],
                GeometricLoading {
                    support_motion_per_factor_m: Some(&[[0.0; 2], [0.1, 0.0]]),
                    ..GeometricLoading::default()
                },
                scales,
                limits(),
                cx
            )
            .is_err()
        );
        let removed = GeometricTruss::try_new(
            &g,
            &case,
            &[0.0],
            &[1000.0],
            GeometricLoading::default(),
            scales,
            limits(),
            cx,
        )
        .unwrap();
        assert_eq!(
            removed.response(&[-1.0], 0.0, cx).unwrap().member_force_n,
            [0.0]
        );
        assert!(
            removed
                .equilibrate(&[[0.0; 2]; 2], 1.0, NewtonKrylovConfig::default(), 5, cx)
                .is_err()
        );
        let state = truss
            .start_from(&[[0.0; 2]; 2], 0.0, &ArcSettings::default(), cx)
            .unwrap();
        let original = state.u.clone();
        let mut malformed = state.clone();
        malformed.u.push(0.0);
        assert!(
            truss
                .advance(&malformed, &ArcSettings::default(), 0, cx)
                .is_err()
        );
        let _ = malformed.u.pop();
        malformed.tangent.0.push(0.0);
        assert!(
            truss
                .advance(&malformed, &ArcSettings::default(), 0, cx)
                .is_err()
        );
        gate.request();
        assert!(matches!(
            truss.advance(&state, &ArcSettings::default(), 1, cx),
            Err(SolidError::Cancelled)
        ));
        assert_eq!(state.u, original);
    });
}
