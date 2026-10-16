//! IRP M1.3-W1a: independent quadratic-solid G0/G1/G2 probes.
//! Native execution belongs to the central batch; no result is claimed merely
//! from these tests' presence. All constitutive inputs are analytic fixtures.

#[path = "tet10/benchmarks.rs"]
mod benchmarks;
#[path = "tet10/convergence.rs"]
mod convergence;
#[path = "tet10/support.rs"]
mod support;

use fs_feec::highorder::simplex::duffy_quadrature;
use fs_solid::linear3::tet10::{
    MAX_QUADRATURE_POINTS_PER_DIR, QuadraticTetProblem, TET10_EDGES, TetQuadrature,
    promote_p1_to_p2, tet10_shapes,
};
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticError, TetLinearElasticProblem, TetMaterialField,
};
use support::*;

#[test]
fn g0_shape_functions_partition_unity_and_kronecker() {
    let mut seed = 58u64;
    for _ in 0..1000 {
        let mut lambda = [0.0; 4];
        for v in &mut lambda {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            *v = ((seed >> 32) as f64 + 1.0) / (u32::MAX as f64 + 1.0);
        }
        let sum: f64 = lambda.iter().sum();
        lambda = lambda.map(|x| x / sum);
        let (n, d) = tet10_shapes(lambda);
        assert!((n.iter().sum::<f64>() - 1.0).abs() < 1e-14);
        for k in 0..3 {
            assert!(d.iter().map(|g| g[k]).sum::<f64>().abs() < 1e-14);
        }
    }
    let mut points = vec![
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    for [a, b] in TET10_EDGES {
        let mut p = [0.0; 4];
        p[a] = 0.5;
        p[b] = 0.5;
        points.push(p);
    }
    for (j, p) in points.into_iter().enumerate() {
        let (n, _) = tet10_shapes(p);
        for (i, v) in n.iter().enumerate() {
            assert_eq!(*v, if i == j { 1.0 } else { 0.0 });
        }
    }
}

fn factorial(n: usize) -> f64 {
    (1..=n).map(|x| x as f64).product()
}

#[test]
fn g0_duffy_quadrature_exactness_matches_rule() {
    for (n, degree) in [(3, 2), (4, 4)] {
        let rule = duffy_quadrature(n);
        assert!(rule.iter().all(|(_, w)| *w > 0.0));
        for a in 0..=degree {
            for b in 0..=degree - a {
                for c in 0..=degree - a - b {
                    let got: f64 = rule
                        .iter()
                        .map(|(l, w)| {
                            w * l[1].powi(a as i32) * l[2].powi(b as i32) * l[3].powi(c as i32)
                        })
                        .sum();
                    let exact =
                        factorial(a) * factorial(b) * factorial(c) / factorial(a + b + c + 3);
                    assert!(
                        (got - exact).abs() < 1e-14,
                        "n={n}, {a},{b},{c}: {got:e} != {exact:e}"
                    );
                }
            }
        }
    }
    let under: f64 = duffy_quadrature(2)
        .iter()
        .map(|(l, w)| w * l[1] * l[1] * l[1] * l[1])
        .sum();
    assert!(
        (under - 1.0 / 210.0).abs() > 1e-6,
        "underintegration falsifier must fail"
    );
}

#[test]
fn g0_consistent_mass_total_and_first_moments_exact() {
    with_cx(|cx| {
        let mesh = support::Mesh {
            nodes: vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            tets: vec![[0, 1, 2, 3]],
        };
        let (nodes, elements) = promote(&mesh, cx);
        let a = assembly(&nodes, &elements, &material(10.0, 0.25), &[], cx);
        fs_la::factor::cholesky(&a.mass.to_dense(), 3 * nodes.len()).unwrap();
        // Total mass and first moments alone cannot detect an incorrect P2
        // mass. Check every entry independently with barycentric moments:
        // integral prod(lambda_i^p_i) = prod(p_i!)/(3+sum(p_i))! on this tet.
        let polynomials: Vec<Vec<([usize; 4], f64)>> = (0..10)
            .map(|i| {
                if i < 4 {
                    let mut square = [0; 4];
                    square[i] = 2;
                    let mut linear = [0; 4];
                    linear[i] = 1;
                    vec![(square, 2.0), (linear, -1.0)]
                } else {
                    let [u, v] = TET10_EDGES[i - 4];
                    let mut product = [0; 4];
                    product[u] = 1;
                    product[v] = 1;
                    vec![(product, 4.0)]
                }
            })
            .collect();
        for i in 0..10 {
            for j in 0..10 {
                let exact: f64 = polynomials[i]
                    .iter()
                    .flat_map(|(p, c)| {
                        polynomials[j].iter().map(move |(q, d)| {
                            let powers: [usize; 4] = core::array::from_fn(|k| p[k] + q[k]);
                            c * d * powers.iter().map(|&k| factorial(k)).product::<f64>()
                                / factorial(3 + powers.iter().sum::<usize>())
                        })
                    })
                    .sum();
                for component in 0..3 {
                    let row = 3 * elements[0][i] + component;
                    let col = 3 * elements[0][j] + component;
                    let (cols, vals) = a.mass.row(row);
                    let got = cols.binary_search(&col).map(|k| vals[k]).unwrap_or(0.0);
                    assert!(
                        (got - exact).abs() < 1e-14,
                        "P2 mass entry ({i},{j},{component}): {got:e} != {exact:e}"
                    );
                }
            }
        }
        for c in 0..3 {
            let one: Vec<_> = (0..3 * nodes.len())
                .map(|i| if i % 3 == c { 1.0 } else { 0.0 })
                .collect();
            let mut mone = vec![0.0; one.len()];
            a.mass.spmv(&one, &mut mone);
            assert!(
                (one.iter().zip(&mone).map(|(x, y)| x * y).sum::<f64>() - 1.0 / 6.0).abs() < 1e-14
            );
            for k in 0..3 {
                let moment: f64 = nodes
                    .iter()
                    .enumerate()
                    .map(|(i, p)| p[k] * mone[3 * i + c])
                    .sum();
                assert!((moment - 1.0 / 24.0).abs() < 1e-14);
            }
        }
        assert!((a.total_mass_kg - 1.0 / 6.0).abs() < 1e-14);
    });
}

#[test]
fn g0_free_patch_has_exactly_six_rigid_modes() {
    with_cx(|cx| {
        let mut mesh = block([2, 1, 1], [1.0; 3]);
        mesh.nodes[4][1] += 0.07;
        mesh.nodes[5][2] -= 0.03;
        let (nodes, elements) = promote(&mesh, cx);
        let a = assembly(&nodes, &elements, &material(10.0, 0.25), &[], cx);
        let modes = fs_modal::eigh_gen_dense(
            &a.stiffness.to_dense(),
            &a.mass.to_dense(),
            a.free_dofs.len(),
        )
        .unwrap();
        let max = modes.last().unwrap().lambda;
        assert_eq!(
            modes
                .iter()
                .filter(|m| m.lambda.abs() <= 1e-10 * max)
                .count(),
            6
        );
        assert!(modes[6].lambda > 1e-6 * max);
    });
}

fn bending(p: [f64; 3], nu: f64) -> [f64; 3] {
    let [x, y, z] = p;
    [-x * y, (x * x + nu * (y * y - z * z)) / 2.0, nu * y * z]
}

#[test]
fn g0_quadratic_patch_is_exact_on_distorted_mesh() {
    with_cx(|cx| {
        let mesh = jittered_cube(3);
        let (nodes, elements) = promote(&mesh, cx);
        let m = material(10.0, 0.25);
        let a = assembly(&nodes, &elements, &m, &[], cx);
        let exact: Vec<_> = nodes.iter().flat_map(|&p| bending(p, 0.25)).collect();
        let boundary: Vec<_> = nodes
            .iter()
            .flat_map(|p| {
                let on = p
                    .iter()
                    .any(|x| x.abs() < 1e-13 || (*x - 1.0).abs() < 1e-13);
                bending(*p, 0.25).map(|u| if on { Some(u) } else { None })
            })
            .collect();
        let u = solve(&a.stiffness, &vec![0.0; exact.len()], &boundary);
        let max = exact.iter().map(|x| x.abs()).fold(0.0, f64::max);
        let error = u
            .iter()
            .zip(&exact)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max)
            / max;
        assert!(
            error < 1e-10,
            "quadratic patch relative displacement error {error:e}"
        );
        let energy = 0.5 * quadratic_form(&a.stiffness, &u);
        assert!(
            (energy / (10.0 / 6.0) - 1.0).abs() < 1e-12,
            "energy {energy:e}"
        );
        let p1 = TetLinearElasticProblem {
            nodes_m: &mesh.nodes,
            tetrahedra: &mesh.tets,
            materials: TetMaterialField::Uniform(&m),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
        }
        .assemble(cx)
        .unwrap();
        let prescribed: Vec<_> = mesh
            .nodes
            .iter()
            .flat_map(|p| {
                let on = p.iter().any(|x| *x == 0.0 || *x == 1.0);
                bending(*p, 0.25).map(|x| if on { Some(x) } else { None })
            })
            .collect();
        let u1 = solve(&p1.stiffness, &vec![0.0; 3 * mesh.nodes.len()], &prescribed);
        let error1 = u1
            .iter()
            .zip(mesh.nodes.iter().flat_map(|&p| bending(p, 0.25)))
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max)
            / max;
        assert!(error1 > 1e-3, "P1 falsifier has no separation: {error1:e}");
    });
}

#[test]
fn g0_rotated_orthotropic_energy_invariant() {
    with_cx(|cx| {
        let mesh = jittered_cube(2);
        let (nodes, elements) = promote(&mesh, cx);
        let law = fs_material::OrthotropicElastic::new(
            [12.0, 5.0, 3.0],
            [0.2, 0.15, 0.1],
            [2.0, 1.5, 1.0],
            1.0,
        )
        .unwrap();
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let rotation = [[0.8, -0.6, 0.0], [0.6, 0.8, 0.0], [0.0, 0.0, 1.0]];
        let rotate =
            |p: [f64; 3]| rotation.map(|r| r.iter().zip(p).map(|(a, b)| a * b).sum::<f64>());
        let m1 = fs_solid::linear3::TetElasticMaterial::try_new_oriented_orthotropic(
            1.0, &law, identity, [1; 32], [2; 32],
        )
        .unwrap();
        let m2 = fs_solid::linear3::TetElasticMaterial::try_new_oriented_orthotropic(
            1.0, &law, rotation, [1; 32], [3; 32],
        )
        .unwrap();
        let a = assembly(&nodes, &elements, &m1, &[], cx);
        let rotated: Vec<_> = nodes.iter().copied().map(rotate).collect();
        let b = assembly(&rotated, &elements, &m2, &[], cx);
        let u: Vec<_> = nodes.iter().flat_map(|&p| bending(p, 0.25)).collect();
        let v: Vec<_> = nodes
            .iter()
            .flat_map(|&p| rotate(bending(p, 0.25)))
            .collect();
        let ea = quadratic_form(&a.stiffness, &u);
        let eb = quadratic_form(&b.stiffness, &v);
        assert!((ea - eb).abs() / ea < 1e-13, "rotation {ea:e} {eb:e}");
    });
}

#[test]
fn g0_refuses_nonpositive_jacobian_at_quadrature_point() {
    with_cx(|cx| {
        let mesh = block([1; 3], [1.0; 3]);
        let (mut nodes, elements) = promote(&mesh, cx);
        nodes[elements[0][4]] = [10.0, -10.0, 10.0];
        let m = material(10.0, 0.25);
        let result = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&m),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
            quadrature: Q,
        }
        .assemble(cx);
        assert!(
            matches!(result, Err(TetElasticError::NonPositiveJacobian { .. })),
            "{result:?}"
        );
    });
}

#[test]
fn g0_refuses_midside_mismatch() {
    with_cx(|cx| {
        let mesh = block([1; 3], [1.0; 3]);
        let (mut nodes, mut elements) = promote(&mesh, cx);
        let mid = elements[0][6];
        nodes.push(nodes[mid]);
        elements[0][6] = nodes.len() - 1;
        let m = material(10.0, 0.25);
        let result = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::Uniform(&m),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
            quadrature: Q,
        }
        .assemble(cx);
        assert!(
            matches!(result, Err(TetElasticError::MidsideMismatch { .. })),
            "{result:?}"
        );
    });
}

#[test]
fn g0_budget_cap_and_cap_plus_one() {
    with_cx(|cx| {
        let mesh = block([1; 3], [1.0; 3]);
        let (nodes, elements) = promote(&mesh, cx);
        let m = material(10.0, 0.25);
        let at = TetAssemblyBudget {
            maximum_nodes: nodes.len(),
            maximum_tetrahedra: elements.len(),
            maximum_free_dofs: 3 * nodes.len(),
            minimum_scaled_jacobian: 1e-12,
        };
        let run = |budget| {
            QuadraticTetProblem {
                nodes_m: &nodes,
                elements: &elements,
                materials: TetMaterialField::Uniform(&m),
                fixed_dofs: &[],
                budget,
                quadrature: Q,
            }
            .assemble(cx)
        };
        run(at).unwrap();
        for under in [
            TetAssemblyBudget {
                maximum_nodes: at.maximum_nodes - 1,
                ..at
            },
            TetAssemblyBudget {
                maximum_tetrahedra: at.maximum_tetrahedra - 1,
                ..at
            },
            TetAssemblyBudget {
                maximum_free_dofs: at.maximum_free_dofs - 1,
                ..at
            },
        ] {
            assert!(
                matches!(run(under),Err(TetElasticError::BudgetExceeded{requested,maximum,..}) if requested==maximum+1)
            );
        }
        TetQuadrature {
            stiffness_points_per_dir: MAX_QUADRATURE_POINTS_PER_DIR,
            mass_points_per_dir: MAX_QUADRATURE_POINTS_PER_DIR,
        }
        .validate()
        .unwrap();
        for q in [
            TetQuadrature {
                stiffness_points_per_dir: 2,
                ..Q
            },
            TetQuadrature {
                mass_points_per_dir: 3,
                ..Q
            },
            TetQuadrature {
                stiffness_points_per_dir: MAX_QUADRATURE_POINTS_PER_DIR + 1,
                ..Q
            },
            TetQuadrature {
                mass_points_per_dir: MAX_QUADRATURE_POINTS_PER_DIR + 1,
                ..Q
            },
        ] {
            assert!(matches!(
                q.validate(),
                Err(TetElasticError::QuadratureInadmissible { .. })
            ));
        }
    });
}

#[test]
fn g0_malformed_inputs_refuse_without_partial_assembly() {
    with_cx(|cx| {
        let mesh = block([1; 3], [1.0; 3]);
        let (nodes, elements) = promote(&mesh, cx);
        let m = material(10.0, 0.25);
        let run = |nodes: &[[f64; 3]], elements: &[[usize; 10]], fixed: &[usize]| {
            QuadraticTetProblem {
                nodes_m: nodes,
                elements,
                materials: TetMaterialField::Uniform(&m),
                fixed_dofs: fixed,
                budget: TetAssemblyBudget::standard(),
                quadrature: Q,
            }
            .assemble(cx)
        };
        assert!(matches!(
            run(&[], &elements, &[]),
            Err(TetElasticError::EmptyMesh)
        ));
        assert!(matches!(
            run(&nodes, &[], &[]),
            Err(TetElasticError::EmptyMesh)
        ));
        let mut bad = nodes.clone();
        bad[0][0] = f64::NAN;
        assert!(matches!(
            run(&bad, &elements, &[]),
            Err(TetElasticError::NonFiniteNode { node: 0 })
        ));
        let mut bad = elements.clone();
        bad[0][0] = nodes.len();
        assert!(matches!(
            run(&nodes, &bad, &[]),
            Err(TetElasticError::NodeIndexOutOfRange { .. })
        ));
        let mut bad = elements.clone();
        bad[0][1] = bad[0][0];
        assert!(matches!(
            run(&nodes, &bad, &[]),
            Err(TetElasticError::RepeatedTetNode { .. })
        ));
        let mut bad = nodes.clone();
        bad.push([2.0; 3]);
        assert!(matches!(
            run(&bad, &elements, &[]),
            Err(TetElasticError::UnreferencedNode { .. })
        ));
        assert!(matches!(
            run(&nodes, &elements, &[0, 0]),
            Err(TetElasticError::DuplicateFixedDof { .. })
        ));
        assert!(matches!(
            run(&nodes, &elements, &[3 * nodes.len()]),
            Err(TetElasticError::FixedDofOutOfRange { .. })
        ));
        assert!(matches!(
            run(&nodes, &elements, &(0..3 * nodes.len()).collect::<Vec<_>>()),
            Err(TetElasticError::NoFreeDofs)
        ));
        assert!(matches!(
            QuadraticTetProblem {
                nodes_m: &nodes,
                elements: &elements,
                materials: TetMaterialField::PerElement(&[]),
                fixed_dofs: &[],
                budget: TetAssemblyBudget::standard(),
                quadrature: Q,
            }
            .assemble(cx),
            Err(TetElasticError::MaterialCountMismatch { .. })
        ));
        // Per-element assignment is the same executable material seam as P1.
        let expanded = vec![m.clone(); elements.len()];
        let a = assembly(&nodes, &elements, &m, &[], cx);
        let b = QuadraticTetProblem {
            nodes_m: &nodes,
            elements: &elements,
            materials: TetMaterialField::PerElement(&expanded),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
            quadrature: Q,
        }
        .assemble(cx)
        .unwrap();
        assert_eq!(a.stiffness, b.stiffness);
        assert_eq!(a.mass, b.mass);
    });
}

#[test]
fn g3_promotion_is_element_permutation_invariant_and_cancellable() {
    with_cx(|cx| {
        let mesh = block([1; 3], [1.0; 3]);
        let (nodes, elements) = promote(&mesh, cx);
        let mut permuted = mesh.clone();
        permuted.tets.reverse();
        let (nodes2, mut elements2) = promote(&permuted, cx);
        elements2.reverse();
        assert_eq!(nodes, nodes2);
        assert_eq!(elements, elements2);
        assert_eq!(
            assembly(&nodes, &elements, &material(10.0, 0.25), &[], cx).stiffness,
            assembly(&nodes, &elements, &material(10.0, 0.25), &[], cx).stiffness
        );
    });
    let gate = fs_exec::CancelGate::new_clock_free();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 58,
                kernel_id: 58,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let mesh = block([1; 3], [1.0; 3]);
        assert!(matches!(
            promote_p1_to_p2(&mesh.nodes, &mesh.tets, None, &cx),
            Err(TetElasticError::Cancelled)
        ));
        let m = material(10.0, 0.25);
        assert!(matches!(
            QuadraticTetProblem {
                nodes_m: &[],
                elements: &[],
                materials: TetMaterialField::Uniform(&m),
                fixed_dofs: &[],
                budget: TetAssemblyBudget::standard(),
                quadrature: Q
            }
            .assemble(&cx),
            Err(TetElasticError::Cancelled)
        ));
    });
}

#[test]
fn g1_mms_orders_isotropic_and_rotated_orthotropic() {
    convergence::mms_orders();
}
#[test]
fn g1_curved_boundary_projection_reduces_geometry_error() {
    convergence::curved_geometry();
}
#[test]
fn g2_cooks_membrane_plane_strain_equivalent() {
    benchmarks::cooks_membrane();
}
#[test]
fn g2_macneal_harder_straight_cantilever() {
    benchmarks::macneal_harder();
}
#[test]
fn g1_free_free_bar_eigen_ladder() {
    benchmarks::free_bar();
}
