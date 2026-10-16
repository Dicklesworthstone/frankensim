//! G0/G1/G2/G3 quadratic hexahedral solids; all geometries are analytic fixtures.
#[path = "hex/benchmarks.rs"]
mod benchmarks;
#[path = "hex/support.rs"]
mod support;
#[path = "tet10/support.rs"]
mod tet_support;
use support::*;

#[test]
fn g0_shape_functions_partition_unity_and_kronecker() {
    for family in FAMILIES {
        let local: &[[i8; 3]] = if family == HexFamily::Serendipity20 {
            &HEX20_REFERENCE_NODES
        } else {
            &HEX27_REFERENCE_NODES
        };
        for (i, p) in local.iter().enumerate() {
            let r = p.map(f64::from);
            let values = if family == HexFamily::Serendipity20 {
                hex20_shapes(r).0.to_vec()
            } else {
                hex27_shapes(r).0.to_vec()
            };
            for (j, v) in values.iter().enumerate() {
                assert!((*v - if i == j { 1.0 } else { 0.0 }).abs() < 1e-14);
            }
        }
        for r in [[0.0; 3], [-0.83, 0.17, 0.51], [0.91, -0.42, -0.67]] {
            let (n, d) = if family == HexFamily::Serendipity20 {
                let (n, d) = hex20_shapes(r);
                (n.to_vec(), d.to_vec())
            } else {
                let (n, d) = hex27_shapes(r);
                (n.to_vec(), d.to_vec())
            };
            assert!((n.iter().sum::<f64>() - 1.0).abs() < 1e-14);
            for k in 0..3 {
                assert!(d.iter().map(|g| g[k]).sum::<f64>().abs() < 1e-14);
            }
        }
    }
}
#[test]
fn g0_gmsh_ordering_permutation() {
    // Independent topology in Gmsh's published local edge and face order.
    let edges = [
        [0, 1],
        [0, 3],
        [0, 4],
        [1, 2],
        [1, 5],
        [2, 3],
        [2, 6],
        [3, 7],
        [4, 5],
        [4, 7],
        [5, 6],
        [6, 7],
    ];
    let faces = [
        [0, 3, 2, 1],
        [0, 1, 5, 4],
        [0, 4, 7, 3],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [4, 5, 6, 7],
    ];
    for (i, [a, b]) in edges.iter().enumerate() {
        let expected = core::array::from_fn::<_, 3, _>(|k| {
            (HEX27_REFERENCE_NODES[*a][k] + HEX27_REFERENCE_NODES[*b][k]) / 2
        });
        assert_eq!(HEX27_REFERENCE_NODES[GMSH27_TO_CANONICAL[8 + i]], expected);
        assert_eq!(GMSH20_TO_CANONICAL[8 + i], GMSH27_TO_CANONICAL[8 + i]);
    }
    for (i, face) in faces.iter().enumerate() {
        let expected = core::array::from_fn::<_, 3, _>(|k| {
            face.iter()
                .map(|v| HEX27_REFERENCE_NODES[*v][k])
                .sum::<i8>()
                / 4
        });
        assert_eq!(HEX27_REFERENCE_NODES[GMSH27_TO_CANONICAL[20 + i]], expected);
    }
    let mut permutation = GMSH27_TO_CANONICAL;
    permutation.sort_unstable();
    assert_eq!(permutation, core::array::from_fn(|i| i));
}
#[test]
fn g0_free_patch_has_exactly_six_rigid_modes() {
    with_cx(|cx| {
        for family in FAMILIES {
            let mesh = mesh([1; 3], [1.0; 3], family, cx);
            let a = assembly(&mesh, &material(10.0, 0.3), &[], cx);
            let modes = fs_modal::eigh_gen_dense(
                &a.stiffness.to_dense(),
                &a.mass.to_dense(),
                a.free_dofs.len(),
            )
            .unwrap();
            let scale = modes.last().unwrap().lambda;
            assert_eq!(
                modes
                    .iter()
                    .filter(|m| m.lambda.abs() < 1e-10 * scale)
                    .count(),
                6
            );
            assert!(modes[6].lambda > 1e-6 * scale);
        }
    });
}
#[test]
fn g0_reduced_integration_refused_and_hourglass_counted() {
    with_cx(|cx| {
        let mesh = mesh([1; 3], [1.0; 3], HexFamily::Serendipity20, cx);
        let law = material(10.0, 0.3);
        let mut p = problem(&mesh, &law, &[]);
        p.quadrature.points_per_dir = 2;
        assert!(matches!(
            p.assemble(cx),
            Err(TetElasticError::UnderIntegrated { points_per_dir: 2 })
        ));
        let HexElements::Hex20(e) = &mesh.elements else {
            unreachable!()
        };
        let nodes = e[0].map(|i| mesh.nodes_m[i]);
        let k = diagnostic_hex20_stiffness(nodes, 2, &law);
        let (values, _) = fs_la::eigen::jacobi_eigh(&k, 60);
        let scale = values.iter().copied().fold(0.0, f64::max);
        let zeros = values.iter().filter(|v| v.abs() < 1e-10 * scale).count();
        assert!(zeros > 6, "reduced HEX20 must show spurious modes");
        println!(
            "{{\"suite\":\"irp-s1\",\"stage\":\"reduced-hex20-diagnostic\",\"zero_modes\":{zeros},\"hourglass_modes\":{}}}",
            zeros - 6
        );
        let full = assembly(&mesh, &law, &[], cx);
        let dense = full.stiffness.to_dense();
        let mut permuted = vec![0.0; 3600];
        for i in 0..60 {
            for j in 0..60 {
                permuted[i * 60 + j] =
                    dense[(3 * e[0][i / 3] + i % 3) * 60 + 3 * e[0][j / 3] + j % 3];
            }
        }
        let independent = diagnostic_hex20_stiffness(nodes, 3, &law);
        let scale = independent.iter().map(|v| v.abs()).fold(0.0, f64::max);
        assert!(
            permuted
                .iter()
                .zip(&independent)
                .all(|(a, b)| (a - b).abs() < 1e-13 * scale)
        );
    });
}
#[test]
fn g0_patch_quadratic_completeness() {
    with_cx(|cx| {
        for family in FAMILIES {
            for distorted in [false, true] {
                let mut mesh = mesh([1; 3], [1.0; 3], family, cx);
                if distorted {
                    for p in &mut mesh.nodes_m {
                        *p = distortion(*p);
                    }
                }
                let mut error = 0.0f64;
                for r in [[0.13, -0.27, 0.41], [-0.45, 0.31, -0.62], [0.7, 0.2, -0.1]] {
                    let (x, n, d, ids) = point_data(&mesh, r);
                    let mut u = [0.0; 3];
                    let mut g = [[0.0; 3]; 3];
                    for a in 0..ids.len() {
                        let v = bending(mesh.nodes_m[ids[a]], 0.3);
                        for i in 0..3 {
                            u[i] += n[a] * v[i];
                            for j in 0..3 {
                                g[i][j] += d[a][j] * v[i];
                            }
                        }
                    }
                    let exact = bending(x, 0.3);
                    let gradient = bending_gradient(x, 0.3);
                    for i in 0..3 {
                        error = error.max((u[i] - exact[i]).abs());
                        for j in 0..3 {
                            error = error.max((g[i][j] - gradient[i][j]).abs());
                        }
                    }
                }
                if !distorted || family == HexFamily::Lagrange27 {
                    assert!(
                        error < 1e-10,
                        "{family:?}, distorted={distorted}, error={error}"
                    );
                }
                println!(
                    "{{\"suite\":\"irp-s1\",\"stage\":\"quadratic-completeness\",\"family\":\"{}\",\"distorted\":{distorted},\"error\":{error},\"gated\":{}}}",
                    family_name(family),
                    !distorted || family == HexFamily::Lagrange27
                );
            }
            for distorted in [false, true] {
                if distorted && family == HexFamily::Serendipity20 {
                    continue;
                }
                let mut mesh = mesh([2; 3], [1.0; 3], family, cx);
                let boundary = mesh
                    .nodes_m
                    .iter()
                    .map(|p| p.iter().any(|v| *v == 0.0 || *v == 1.0))
                    .collect::<Vec<_>>();
                if distorted {
                    for p in &mut mesh.nodes_m {
                        *p = distortion(*p);
                    }
                }
                let a = assembly(&mesh, &material(10.0, 0.3), &[], cx);
                let prescribed = mesh
                    .nodes_m
                    .iter()
                    .zip(&boundary)
                    .flat_map(|(p, boundary)| {
                        bending(*p, 0.3).map(|v| if *boundary { Some(v) } else { None })
                    })
                    .collect::<Vec<_>>();
                let u = super_solve(&a, &prescribed);
                let exact = mesh
                    .nodes_m
                    .iter()
                    .flat_map(|p| bending(*p, 0.3))
                    .collect::<Vec<_>>();
                assert!(u.iter().zip(exact).all(|(a, b)| (a - b).abs() < 1e-10));
            }
        }
    });
}
fn super_solve(a: &TetElasticAssembly, prescribed: &[Option<f64>]) -> Vec<f64> {
    tet_support::solve(&a.stiffness, &vec![0.0; prescribed.len()], prescribed)
}
#[test]
fn g0_consistent_mass_total_and_first_moments_exact() {
    with_cx(|cx| {
        for family in FAMILIES {
            let mesh = mesh([1; 3], [1.0; 3], family, cx);
            let a = assembly(&mesh, &material(10.0, 0.3), &[], cx);
            let size = a.free_dofs.len();
            fs_la::factor::cholesky(&a.mass.to_dense(), size).unwrap();
            assert!((a.total_mass_kg - 1.0).abs() < 1e-14);
            for c in 0..3 {
                let one = (0..size)
                    .map(|i| if i % 3 == c { 1.0 } else { 0.0 })
                    .collect::<Vec<_>>();
                let mut mass = vec![0.0; size];
                a.mass.spmv(&one, &mut mass);
                if let HexElements::Hex20(e) = &mesh.elements {
                    assert!(
                        mass[3 * e[0][0] + c] < 0.0,
                        "HEX20 row-sum lumping has negative corner masses"
                    );
                }
                assert!(
                    (one.iter().zip(&mass).map(|(a, b)| a * b).sum::<f64>() - 1.0).abs() < 1e-14
                );
                for k in 0..3 {
                    let moment = mesh
                        .nodes_m
                        .iter()
                        .enumerate()
                        .map(|(i, p)| p[k] * mass[3 * i + c])
                        .sum::<f64>();
                    assert!((moment - 0.5).abs() < 1e-14);
                }
            }
            let law = material(10.0, 0.3);
            let mut higher = problem(&mesh, &law, &[]);
            higher.quadrature.points_per_dir = 4;
            let high = higher.assemble(cx).unwrap();
            let high = high.mass.to_dense();
            let low = a.mass.to_dense();
            let scale = high.iter().map(|v| v.abs()).fold(0.0, f64::max);
            assert!(
                low.iter()
                    .zip(&high)
                    .all(|(a, b)| (a - b).abs() <= 1e-13 * scale)
            );
        }
    });
}
#[test]
fn g0_refuses_nonpositive_jacobian_at_gauss_point() {
    with_cx(|cx| {
        for family in FAMILIES {
            let mut mesh = mesh([1; 3], [1.0; 3], family, cx);
            for p in &mut mesh.nodes_m {
                p[0] = -p[0];
            }
            assert!(matches!(
                problem(&mesh, &material(10.0, 0.3), &[]).assemble(cx),
                Err(TetElasticError::NonPositiveJacobian { .. })
            ));
            let mut folded = mesh;
            for p in &mut folded.nodes_m {
                // Undo the preceding reflection, then create the trilinear
                // map (r,s,(1-2r)t). Its centre and first 18 Gauss points are
                // positive; the final x slab is inverted.
                let r = -2.0 * p[0] - 1.0;
                let s = 2.0 * p[1] - 1.0;
                let t = 2.0 * p[2] - 1.0;
                *p = [r, s, (1.0 - 2.0 * r) * t];
            }
            assert!(matches!(
                problem(&folded, &material(10.0, 0.3), &[]).assemble(cx),
                Err(TetElasticError::NonPositiveJacobian { point: 18, .. })
            ));
        }
    });
}
#[test]
fn g0_budget_cap_and_cap_plus_one() {
    assert!(HexQuadrature { points_per_dir: 12 }.validate().is_ok());
    assert!(matches!(
        HexQuadrature { points_per_dir: 13 }.validate(),
        Err(TetElasticError::QuadratureInadmissible { .. })
    ));
    with_cx(|cx| {
        for family in FAMILIES {
            let mesh = mesh([1; 3], [1.0; 3], family, cx);
            let law = material(10.0, 0.3);
            let mut p = problem(&mesh, &law, &[]);
            let dofs = 3 * mesh.nodes_m.len();
            p.budget.maximum_free_dofs = dofs;
            p.budget.maximum_tetrahedra = 1;
            assert!(p.assemble(cx).is_ok());
            p.budget.maximum_free_dofs = dofs - 1;
            assert!(matches!(
                p.assemble(cx),
                Err(TetElasticError::BudgetExceeded {
                    what: "free_dofs",
                    ..
                })
            ));
            p.budget.maximum_free_dofs = dofs;
            p.budget.maximum_tetrahedra = 0;
            assert!(matches!(
                p.assemble(cx),
                Err(TetElasticError::BudgetExceeded {
                    what: "tetrahedra",
                    ..
                })
            ));
        }
    });
}
#[test]
fn g0_refuses_mixed_family_connectivity() {
    with_cx(|cx| {
        for family in FAMILIES {
            let mesh = mesh([1; 3], [1.0; 3], family, cx);
            let law = material(10.0, 0.3);
            let mut p = problem(&mesh, &law, &[]);
            p.family = if family == HexFamily::Serendipity20 {
                HexFamily::Lagrange27
            } else {
                HexFamily::Serendipity20
            };
            assert!(matches!(
                p.assemble(cx),
                Err(TetElasticError::HexConnectivityMismatch { .. })
            ));
        }
    });
}

#[test]
fn g4_cancellation_refuses_before_assembly_or_mesh_publication() {
    let mesh = with_cx(|cx| mesh([1; 3], [1.0; 3], HexFamily::Lagrange27, cx));
    let gate = fs_exec::CancelGate::new_clock_free();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 59,
                kernel_id: 59,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        assert!(matches!(
            problem(&mesh, &material(1.0, 0.3), &[]).assemble(&cx),
            Err(TetElasticError::Cancelled)
        ));
        assert!(matches!(
            extruded_box_hex(
                [1; 3],
                [1.0; 3],
                HexFamily::Lagrange27,
                HexMeshBudget {
                    maximum_nodes: 27,
                    maximum_elements: 1
                },
                &cx
            ),
            Err(fs_mesh::structured_hex::StructuredHexError::Cancelled)
        ));
    });
}
#[test]
fn g1_roller_block_exact_p_modes_order() {
    benchmarks::roller();
}
#[test]
fn g2_slender_square_bar_free_free_vs_timoshenko() {
    benchmarks::slender_bar();
}
#[test]
fn g3_hex27_vs_tet10_same_bar() {
    benchmarks::cross_family();
}
