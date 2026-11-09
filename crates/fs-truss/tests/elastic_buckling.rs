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

// G1/G3 total critical-load derivatives through both compatible physical states.
#[allow(clippy::too_many_arguments)]
fn design_factor(
    nodes: &[[f64; 2]],
    reference_loads: &[[f64; 2]],
    areas: &[f64],
    young: &[f64],
    strain: &[f64],
    motion: &[[f64; 2]],
    fixed: &[[f64; 2]],
    cx: &Cx<'_>,
) -> f64 {
    let members = [(0, 1), (1, 2), (1, 3)];
    let lengths: Vec<_> = members
        .iter()
        .map(|&(a, b)| (nodes[b][0] - nodes[a][0]).hypot(nodes[b][1] - nodes[a][1]))
        .collect();
    let ground =
        GroundStructure::try_from_parts(nodes, &members, &lengths, GroundLimits::default(), cx)
            .unwrap();
    let case = LayoutCase::try_new(
        vec![[true; 2], [false; 2], [true; 2], [true; 2]],
        reference_loads.to_vec(),
        4,
    )
    .unwrap();
    solve(
        &ground,
        &case,
        areas,
        young,
        BucklingLoadFamily {
            imposed: ElasticImposedState {
                member_stress_free_strain: Some(strain),
                prescribed_displacements_m: Some(motion),
            },
            fixed_loads_n: Some(fixed),
        },
        BucklingControls::default(),
        cx,
    )
    .unwrap()
    .buckling
    .loads[0]
}

fn derivative_near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-5 * actual.abs().max(expected.abs()).max(1e-4),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_buckling_design_derivatives_keep_the_fixed_and_reference_loads_distinct() {
    use fs_truss::elastic::buckling::sensitivity::{
        BucklingSensitivityControls, PreparedBucklingDesign,
    };
    with_cx(|cx, _| {
        let (ground, case) = braced_node(cx, 0.0, 1.0);
        let areas = [1.0; 3];
        let young = [2e6, 1e6, 100.0];
        let strain = [1e-5, 1e-5, 0.0];
        let motion = [[0.0; 2], [0.0; 2], [2e-5, 0.0], [0.0; 2]];
        let fixed = [[0.0; 2], [20.0, 0.0], [0.0; 2], [0.0; 2]];
        let design = PreparedBucklingDesign::try_new(
            &ground,
            &case,
            &areas,
            &young,
            BucklingLoadFamily {
                imposed: ElasticImposedState {
                    member_stress_free_strain: Some(&strain),
                    prescribed_displacements_m: Some(&motion),
                },
                fixed_loads_n: Some(&fixed),
            },
            BucklingControls::default(),
            cx,
        )
        .unwrap();
        assert_eq!(design.retained_modes(), 1);
        let result = design
            .solve_mode(0, BucklingSensitivityControls::default(), cx)
            .unwrap();
        near(result.load_factor, 350.0);
        near(result.geometric_modal_work, 1.0);
        near(result.fixed_load_gradient[1][0], -1.0);
        near(result.reference_load_gradient[1][0], -350.0);
        near(result.prescribed_displacement_gradient[0][0], -3e6);
        near(result.prescribed_displacement_gradient[2][0], 3e6);
        near(
            result.members[0].unwrap().stress_free_strain_derivative,
            -6e6,
        );
        near(
            result.members[1].unwrap().stress_free_strain_derivative,
            -3e6,
        );
        near(result.members[2].unwrap().young_modulus_derivative, 4.0);
        near(result.members[2].unwrap().area_derivative, 400.0);
        assert!(result.relative_eigen_residual <= 1e-8);
        assert!(result.adjoint_relative_equilibrium_residual < 1e-9);
        for node in [0, 2, 3] {
            assert_eq!(result.fixed_load_gradient[node], [0.0; 2]);
            assert_eq!(result.reference_load_gradient[node], [0.0; 2]);
        }
        assert_eq!(result.prescribed_displacement_gradient[1], [0.0; 2]);
        let again = design
            .solve_mode(0, BucklingSensitivityControls::default(), cx)
            .unwrap();
        assert_eq!(result.members, again.members);
        assert_eq!(
            result.node_coordinate_gradient,
            again.node_coordinate_gradient
        );
    });
}

#[test]
#[allow(clippy::too_many_lines)]
fn g3_buckling_design_adjoints_match_member_load_support_and_geometry_differences() {
    use fs_truss::elastic::buckling::sensitivity::{
        BucklingSensitivityControls, PreparedBucklingDesign,
    };
    with_cx(|cx, _| {
        let angle = 0.37_f64;
        let (sn, cs) = angle.sin_cos();
        let rotate = |v: [f64; 2]| [cs * v[0] - sn * v[1], sn * v[0] + cs * v[1]];
        let (ground, case) = braced_node(cx, angle, 1.3);
        let areas = [1.0; 3];
        let young = [2e6, 1e6, 100.0];
        let strain = [1e-5, 2e-5, 0.0];
        let motion = [[-1e-5, 0.0], [0.0; 2], [3e-5, 0.0], [0.0, 1e-5]].map(rotate);
        let fixed = [[0.0; 2], [20.0, 0.0], [0.0; 2], [0.0; 2]].map(rotate);
        let design = PreparedBucklingDesign::try_new(
            &ground,
            &case,
            &areas,
            &young,
            BucklingLoadFamily {
                imposed: ElasticImposedState {
                    member_stress_free_strain: Some(&strain),
                    prescribed_displacements_m: Some(&motion),
                },
                fixed_loads_n: Some(&fixed),
            },
            BucklingControls::default(),
            cx,
        )
        .unwrap();
        let result = design
            .solve_mode(0, BucklingSensitivityControls::default(), cx)
            .unwrap();
        for index in 0..3 {
            for field in 0..3 {
                let (mut ap, mut am, mut ep, mut em, mut sp, mut sm) =
                    (areas, areas, young, young, strain, strain);
                let derivative = result.members[index].unwrap();
                let (h, predicted) = match field {
                    0 => {
                        let h = 1e-5;
                        ap[index] += h;
                        am[index] -= h;
                        (h, derivative.area_derivative)
                    }
                    1 => {
                        let h = young[index] * 1e-5;
                        ep[index] += h;
                        em[index] -= h;
                        (h, derivative.young_modulus_derivative)
                    }
                    _ => {
                        let h = 1e-9;
                        sp[index] += h;
                        sm[index] -= h;
                        (h, derivative.stress_free_strain_derivative)
                    }
                };
                let plus = design_factor(
                    ground.nodes(),
                    case.loads(),
                    &ap,
                    &ep,
                    &sp,
                    &motion,
                    &fixed,
                    cx,
                );
                let minus = design_factor(
                    ground.nodes(),
                    case.loads(),
                    &am,
                    &em,
                    &sm,
                    &motion,
                    &fixed,
                    cx,
                );
                derivative_near(predicted, (plus - minus) / (2.0 * h));
            }
        }
        for node in 0..4 {
            for axis in 0..2 {
                let h = 2e-6;
                let (mut xp, mut xm) = (ground.nodes().to_vec(), ground.nodes().to_vec());
                xp[node][axis] += h;
                xm[node][axis] -= h;
                derivative_near(
                    result.node_coordinate_gradient[node][axis],
                    (design_factor(
                        &xp,
                        case.loads(),
                        &areas,
                        &young,
                        &strain,
                        &motion,
                        &fixed,
                        cx,
                    ) - design_factor(
                        &xm,
                        case.loads(),
                        &areas,
                        &young,
                        &strain,
                        &motion,
                        &fixed,
                        cx,
                    )) / (2.0 * h),
                );
                if node == 1 {
                    let (mut fp, mut fm) = (fixed, fixed);
                    let h = 1e-4;
                    fp[node][axis] += h;
                    fm[node][axis] -= h;
                    derivative_near(
                        result.fixed_load_gradient[node][axis],
                        (design_factor(
                            ground.nodes(),
                            case.loads(),
                            &areas,
                            &young,
                            &strain,
                            &motion,
                            &fp,
                            cx,
                        ) - design_factor(
                            ground.nodes(),
                            case.loads(),
                            &areas,
                            &young,
                            &strain,
                            &motion,
                            &fm,
                            cx,
                        )) / (2.0 * h),
                    );
                    let (mut rp, mut rm) = (case.loads().to_vec(), case.loads().to_vec());
                    let h = 1e-6;
                    rp[node][axis] += h;
                    rm[node][axis] -= h;
                    derivative_near(
                        result.reference_load_gradient[node][axis],
                        (design_factor(
                            ground.nodes(),
                            &rp,
                            &areas,
                            &young,
                            &strain,
                            &motion,
                            &fixed,
                            cx,
                        ) - design_factor(
                            ground.nodes(),
                            &rm,
                            &areas,
                            &young,
                            &strain,
                            &motion,
                            &fixed,
                            cx,
                        )) / (2.0 * h),
                    );
                } else {
                    let (mut up, mut um) = (motion, motion);
                    let h = 1e-9;
                    up[node][axis] += h;
                    um[node][axis] -= h;
                    derivative_near(
                        result.prescribed_displacement_gradient[node][axis],
                        (design_factor(
                            ground.nodes(),
                            case.loads(),
                            &areas,
                            &young,
                            &strain,
                            &up,
                            &fixed,
                            cx,
                        ) - design_factor(
                            ground.nodes(),
                            case.loads(),
                            &areas,
                            &young,
                            &strain,
                            &um,
                            &fixed,
                            cx,
                        )) / (2.0 * h),
                    );
                }
            }
        }
        let mut reference_scale = 0.0;
        let mut shape_scale = 0.0;
        let mut rotation = 0.0;
        for node in 0..4 {
            let turn = |v: [f64; 2]| [-v[1], v[0]];
            for axis in 0..2 {
                reference_scale +=
                    result.reference_load_gradient[node][axis] * case.loads()[node][axis];
                shape_scale += result.node_coordinate_gradient[node][axis]
                    * ground.nodes()[node][axis]
                    + result.prescribed_displacement_gradient[node][axis] * motion[node][axis];
                rotation += result.node_coordinate_gradient[node][axis]
                    * turn(ground.nodes()[node])[axis]
                    + result.prescribed_displacement_gradient[node][axis]
                        * turn(motion[node])[axis]
                    + result.fixed_load_gradient[node][axis] * turn(fixed[node])[axis]
                    + result.reference_load_gradient[node][axis] * turn(case.loads()[node])[axis];
            }
        }
        near(reference_scale, -result.load_factor);
        // No mass: uniformly scaling reference geometry and imposed motion
        // scales both pencil matrices identically, leaving lambda unchanged.
        assert!(shape_scale.abs() < 1e-7 * result.load_factor);
        assert!(rotation.abs() < 1e-7 * result.load_factor);
        for axis in 0..2 {
            let sum: f64 = result
                .node_coordinate_gradient
                .iter()
                .map(|g| g[axis])
                .sum();
            assert!(sum.abs() < 1e-10 * result.load_factor);
        }
    });
}

#[test]
fn g4_buckling_derivatives_refuse_clusters_guards_empty_roots_work_and_cancellation() {
    use fs_truss::elastic::buckling::sensitivity::{
        BucklingSensitivityControls, PreparedBucklingDesign,
    };
    with_cx(|cx, gate| {
        let nodes = [
            [-2.0, 0.0],
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, -1.0],
            [-2.0, 4.0],
            [0.0, 4.0],
            [1.0, 4.0],
            [0.0, 3.0],
        ];
        let ground = GroundStructure::try_from_parts(
            &nodes,
            &[(0, 1), (1, 2), (1, 3), (4, 5), (5, 6), (5, 7)],
            &[2.0, 1.0, 1.0, 2.0, 1.0, 1.0],
            GroundLimits::default(),
            cx,
        )
        .unwrap();
        let supports = vec![
            [true; 2], [false; 2], [true; 2], [true; 2], [true; 2], [false; 2], [true; 2],
            [true; 2],
        ];
        let mut loads = vec![[0.0; 2]; 8];
        loads[1] = [1.0, 0.0];
        loads[5] = [1.0, 0.0];
        let case = LayoutCase::try_new(supports.clone(), loads.clone(), 8).unwrap();
        let areas = [1.0; 6];
        let young = [2e6, 1e6, 100.0, 2e6, 1e6, 100.0];
        let cluster = PreparedBucklingDesign::try_new(
            &ground,
            &case,
            &areas,
            &young,
            BucklingLoadFamily::default(),
            BucklingControls {
                modes: 2,
                ..Default::default()
            },
            cx,
        )
        .unwrap();
        assert_eq!(cluster.retained_modes(), 2);
        for index in 0..2 {
            assert!(
                cluster
                    .solve_mode(index, BucklingSensitivityControls::default(), cx)
                    .is_err()
            );
        }
        loads[5] = [0.7, 0.0];
        let case = LayoutCase::try_new(supports.clone(), loads, 8).unwrap();
        let design = PreparedBucklingDesign::try_new(
            &ground,
            &case,
            &areas,
            &young,
            BucklingLoadFamily::default(),
            BucklingControls::default(),
            cx,
        )
        .unwrap();
        assert_eq!(design.retained_modes(), 1);
        assert!(design.analysis().buckling.loads.len() >= 2);
        let accepted = design
            .solve_mode(0, BucklingSensitivityControls::default(), cx)
            .unwrap();
        assert!(accepted.returned_root_relative_gap.unwrap() > 0.1);
        assert!(
            design
                .solve_mode(1, BucklingSensitivityControls::default(), cx)
                .is_err()
        );
        assert!(
            design
                .solve_mode(
                    0,
                    BucklingSensitivityControls {
                        maximum_bytes: 1,
                        ..Default::default()
                    },
                    cx
                )
                .is_err()
        );
        assert!(
            design
                .solve_mode(
                    0,
                    BucklingSensitivityControls {
                        elastic: ElasticControls {
                            maximum_iterations: 0,
                            ..Default::default()
                        },
                        ..Default::default()
                    },
                    cx
                )
                .is_err()
        );
        let empty_case = LayoutCase::try_new(supports, vec![[0.0; 2]; 8], 8).unwrap();
        let empty = PreparedBucklingDesign::try_new(
            &ground,
            &empty_case,
            &areas,
            &young,
            BucklingLoadFamily::default(),
            BucklingControls::default(),
            cx,
        )
        .unwrap();
        assert_eq!(empty.retained_modes(), 0);
        assert!(
            empty
                .solve_mode(0, BucklingSensitivityControls::default(), cx)
                .is_err()
        );
        gate.request();
        assert!(
            design
                .solve_mode(0, BucklingSensitivityControls::default(), cx)
                .is_err()
        );
    });
}
