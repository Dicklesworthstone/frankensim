//! G0/G1/G3/G4 native HEX8 tensor stiffness, consistent inertia and consumers.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::dynamics::{ExactLinearSolid, ExactLinearSolidControls};
use fs_solid::linear::Formulation;
use fs_solid::linear3::hex::hex8_point;
use fs_solid::linear3::hex8::Hex8ElasticProblem;
use fs_solid::linear3::{TetAssemblyBudget, TetElasticError, TetElasticMaterial, TetMaterialField};

const CORNERS: [[f64; 3]; 8] = [
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
];
const CELLS: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];
const FORMS: [Formulation; 2] = [Formulation::Standard, Formulation::BBar];

fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 81,
                kernel_id: 81,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn coupled_material() -> TetElasticMaterial {
    let diagonal = [11.0, 14.0, 9.0, 6.0, 5.0, 8.0];
    let a = [0.8, -0.4, 0.3, 0.7, -0.6, 0.2];
    let b = [-0.2, 0.5, 0.9, -0.3, 0.4, 0.6];
    let c = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (if i == j { diagonal[i] } else { 0.0 }) + a[i] * a[j] + b[i] * b[j]
        })
    });
    TetElasticMaterial::try_new_mandel(2.5, c, [81; 32]).unwrap()
}

fn problem<'a>(
    nodes: &'a [[f64; 3]],
    material: &'a TetElasticMaterial,
    formulation: Formulation,
) -> Hex8ElasticProblem<'a> {
    Hex8ElasticProblem {
        nodes_m: nodes,
        hexahedra: &CELLS,
        materials: TetMaterialField::Uniform(material),
        formulation,
        fixed_dofs: &[],
        budget: TetAssemblyBudget::standard(),
        stiffness_points_per_dir: 2,
    }
}

fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "actual={actual:e}, expected={expected:e}, allowance={allowance:e}"
    );
}

fn work(matrix: &[f64], left: &[f64], right: &[f64]) -> f64 {
    let n = left.len();
    (0..n)
        .map(|i| left[i] * (0..n).map(|j| matrix[i * n + j] * right[j]).sum::<f64>())
        .sum()
}

fn constitutive_work(material: &TetElasticMaterial, strain: [f64; 6]) -> f64 {
    let c = material.stiffness_mandel_pa();
    (0..6)
        .map(|i| strain[i] * (0..6).map(|j| c[i][j] * strain[j]).sum::<f64>())
        .sum()
}

fn warped() -> [[f64; 3]; 8] {
    // Exact trilinear map: x=r+a*s*t, y=s+b*r*t, z=t+c*r*s.
    CORNERS.map(|[r, s, t]| [r + 0.12 * s * t, s + 0.08 * r * t, t + 0.06 * r * s])
}

#[test]
fn g0_g3_affine_tensor_patch_and_rigid_motion_on_warped_bricks() {
    with_cx(false, |cx| {
        let material = coupled_material();
        let nodes = warped();
        let a = [
            [0.012, -0.008, 0.004],
            [0.006, -0.007, 0.003],
            [-0.005, 0.009, 0.011],
        ];
        let root_two = 2.0_f64.sqrt();
        let strain = [
            a[0][0],
            a[1][1],
            a[2][2],
            (a[0][1] + a[1][0]) / root_two,
            (a[1][2] + a[2][1]) / root_two,
            (a[2][0] + a[0][2]) / root_two,
        ];
        let u: Vec<_> = nodes
            .iter()
            .flat_map(|x| a.map(|row| 0.013 + (0..3).map(|j| row[j] * x[j]).sum::<f64>()))
            .collect();
        // det(J)=1-ab*t^2-ac*s^2-bc*r^2+2abc*r*s*t.
        let volume = 8.0 * (1.0 - (0.12 * 0.08 + 0.12 * 0.06 + 0.08 * 0.06) / 3.0);
        let expected = volume * constitutive_work(&material, strain);
        for form in FORMS {
            for order in [2, 4] {
                let mut p = problem(&nodes, &material, form);
                p.stiffness_points_per_dir = order;
                let assembly = p.assemble(cx).unwrap();
                let k = assembly.stiffness.to_dense();
                near(assembly.element_volumes_m3[0], volume, 3e-13);
                near(work(&k, &u, &u), expected, 3e-12 * expected);
                let scale = k.iter().copied().map(f64::abs).fold(0.0, f64::max);
                for axis in 0..3 {
                    let translation: Vec<_> = (0..24)
                        .map(|d| if d % 3 == axis { 1.0 } else { 0.0 })
                        .collect();
                    let rotation: Vec<_> = nodes
                        .iter()
                        .flat_map(|x| {
                            std::array::from_fn::<_, 3, _>(|i| {
                                if i == (axis + 1) % 3 {
                                    -x[(axis + 2) % 3]
                                } else if i == (axis + 2) % 3 {
                                    x[(axis + 1) % 3]
                                } else {
                                    0.0
                                }
                            })
                        })
                        .collect();
                    for rigid in [translation, rotation] {
                        for row in k.chunks_exact(24) {
                            let residual: f64 = row.iter().zip(&rigid).map(|(x, y)| x * y).sum();
                            near(residual, 0.0, 4e-12 * scale);
                        }
                    }
                }
            }
        }
    });
}

#[test]
fn g1_full_anisotropic_bbar_work_matches_an_independent_polynomial_integral() {
    with_cx(false, |cx| {
        let material = coupled_material();
        let [a, b, c, d, e, f] = [0.013, -0.009, 0.011, 0.007, -0.005, 0.003];
        let u: Vec<_> = CORNERS
            .iter()
            .flat_map(|[x, y, z]| {
                [
                    a * x * y + b * y * z,
                    c * y * z + d * x * z,
                    e * z * x + f * x * y,
                ]
            })
            .collect();
        let s = std::f64::consts::FRAC_1_SQRT_2;
        // These are the exact Mandel coefficients of x, y and z. Their cross
        // integrals vanish on [-1,1]^3 and each square integrates to 8/3.
        let coefficients = [
            [0.0, 0.0, e, a * s, (d + f) * s, 0.0],
            [a, 0.0, 0.0, 0.0, c * s, (f + b) * s],
            [0.0, c, 0.0, (b + d) * s, 0.0, e * s],
        ];
        for form in FORMS {
            let expected: f64 = coefficients
                .iter()
                .map(|&coefficient| {
                    let mut strain = coefficient;
                    if form == Formulation::BBar {
                        // The exact volume-mean dilatation is zero.
                        let mean = (strain[0] + strain[1] + strain[2]) / 3.0;
                        for value in &mut strain[..3] {
                            *value -= mean;
                        }
                    }
                    (8.0 / 3.0) * constitutive_work(&material, strain)
                })
                .sum();
            let assembly = problem(&CORNERS, &material, form).assemble(cx).unwrap();
            near(
                work(&assembly.stiffness.to_dense(), &u, &u),
                expected,
                2e-12 * expected,
            );
        }
    });
}

#[test]
fn g1_consistent_mass_retains_warped_geometry_and_is_independent_of_bbar() {
    with_cx(false, |cx| {
        let material = coupled_material();
        let standard = problem(&CORNERS, &material, Formulation::Standard)
            .assemble(cx)
            .unwrap();
        let bbar = problem(&CORNERS, &material, Formulation::BBar)
            .assemble(cx)
            .unwrap();
        let mass = standard.mass.to_dense();
        assert_eq!(mass, bbar.mass.to_dense());
        near(standard.total_mass_kg, 20.0, 2e-13);
        for a in 0..8 {
            for b in 0..8 {
                let scalar = material.density_kg_m3
                    * (0..3)
                        .map(|i| {
                            if CORNERS[a][i] == CORNERS[b][i] {
                                2.0 / 3.0
                            } else {
                                1.0 / 3.0
                            }
                        })
                        .product::<f64>();
                for i in 0..3 {
                    for j in 0..3 {
                        near(
                            mass[(3 * a + i) * 24 + 3 * b + j],
                            if i == j { scalar } else { 0.0 },
                            2e-14,
                        );
                    }
                }
            }
        }
        let spin: Vec<_> = CORNERS.iter().flat_map(|[x, y, _]| [-y, *x, 0.0]).collect();
        near(
            0.5 * work(&mass, &spin, &spin),
            8.0 * material.density_kg_m3 / 3.0,
            2e-13,
        );

        let nodes = warped();
        let p = problem(&nodes, &material, Formulation::BBar);
        let warped = p.assemble(cx).unwrap();
        let [a, b, c] = [0.12, 0.08, 0.06];
        // Independent exact integral of N_0^2 det(J), using
        // integral (1-r)^2=8/3, integral r(1-r)^2=-4/3,
        // integral r^2(1-r)^2=16/15. A two-point mass rule fails this.
        let exact_diagonal = material.density_kg_m3
            * (8.0 / 27.0 - (a * b + a * c + b * c) * 16.0 / 135.0 - 2.0 * a * b * c / 27.0);
        near(warped.mass.to_dense()[0], exact_diagonal, 3e-14);
        let mut refined = p;
        refined.stiffness_points_per_dir = 5;
        let refined = refined.assemble(cx).unwrap();
        for (actual, expected) in warped.mass.to_dense().iter().zip(refined.mass.to_dense()) {
            near(*actual, expected, 3e-14);
        }
    });
}

#[test]
fn g1_bbar_removes_the_pointwise_bulk_penalty_without_removing_shear() {
    with_cx(false, |cx| {
        let shear = 3.0;
        let amplitude = 0.001;
        let u: Vec<_> = CORNERS
            .iter()
            .flat_map(|[x, y, _]| [amplitude * x * y, 0.0, 0.0])
            .collect();
        // u_x=a*x*y has zero mean divergence. Its retained deviatoric strain
        // gives integral 2G dev(eps):dev(eps) = 56G*a^2/9.
        let expected_bbar = 56.0 * shear * amplitude * amplitude / 9.0;
        for bulk in [3.0, 3.0e6] {
            let c = std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    if i < 3 && j < 3 {
                        bulk + if i == j {
                            4.0 * shear / 3.0
                        } else {
                            -2.0 * shear / 3.0
                        }
                    } else if i == j {
                        2.0 * shear
                    } else {
                        0.0
                    }
                })
            });
            let material = TetElasticMaterial::try_new_mandel(2.5, c, [82; 32]).unwrap();
            for form in FORMS {
                let assembly = problem(&CORNERS, &material, form).assemble(cx).unwrap();
                let actual = work(&assembly.stiffness.to_dense(), &u, &u);
                let expected = if form == Formulation::Standard {
                    expected_bbar + 8.0 * bulk * amplitude * amplitude / 3.0
                } else {
                    expected_bbar
                };
                near(actual, expected, 3e-8 * expected);
            }
        }
    });
}

#[test]
fn g1_native_modal_consumer_retains_six_rigid_modes_and_no_hourglass_modes() {
    with_cx(false, |cx| {
        let nodes = warped();
        let material = coupled_material();
        for form in FORMS {
            let assembly = problem(&nodes, &material, form).assemble(cx).unwrap();
            let modes = fs_modal::eigh_gen_dense(
                &assembly.stiffness.to_dense(),
                &assembly.mass.to_dense(),
                assembly.free_dofs.len(),
            )
            .unwrap();
            assert_eq!(modes.len(), 24);
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
fn g1_native_time_consumer_uses_the_same_supported_mass_and_stiffness() {
    with_cx(false, |cx| {
        let material = coupled_material();
        let fixed: Vec<_> = (0..24).filter(|d| *d != 18).collect();
        let mut p = problem(&CORNERS, &material, Formulation::BBar);
        p.fixed_dofs = &fixed;
        let assembly = p.assemble(cx).unwrap();
        assert_eq!(assembly.free_dofs, vec![18]);
        let stiffness = assembly.stiffness.to_dense()[0];
        let mass = assembly.mass.to_dense()[0];
        near(mass, 20.0 / 27.0, 2e-14);
        let omega = (stiffness / mass).sqrt();
        let phase = 0.2_f64;
        let (u, v) = (0.001, -0.002);
        let mut solid = ExactLinearSolid::prepare(
            &assembly,
            phase / omega,
            &[u],
            &[v],
            ExactLinearSolidControls {
                maximum_free_dofs: 1,
                maximum_preparation_bytes: 64 * 1024 * 1024,
                maximum_exponential_multiply_terms: 100_000_000,
                maximum_relative_operator_error: 1e-12,
                maximum_work_error_j: 1e-12,
                maximum_steps: 1,
            },
            cx,
        )
        .unwrap();
        let energy = 0.5 * (stiffness * u * u + mass * v * v);
        solid.advance(&[0.0], cx).unwrap();
        near(
            solid.displacement_m()[0],
            u * phase.cos() + v / omega * phase.sin(),
            3e-14,
        );
        near(
            solid.velocity_mps()[0],
            -omega * u * phase.sin() + v * phase.cos(),
            3e-14,
        );
        near(
            solid.stored_energy_j() + solid.kinetic_energy_j(),
            energy,
            1e-13,
        );
    });
}

#[test]
fn g0_g4_geometry_constraints_quadrature_budgets_and_cancellation_refuse() {
    with_cx(false, |cx| {
        let material = coupled_material();
        let p = problem(&CORNERS, &material, Formulation::BBar);
        let mut bad = p;
        bad.stiffness_points_per_dir = 1;
        assert!(matches!(
            bad.assemble(cx),
            Err(TetElasticError::UnderIntegrated { .. })
        ));
        bad.stiffness_points_per_dir = 13;
        assert!(matches!(
            bad.assemble(cx),
            Err(TetElasticError::QuadratureInadmissible { .. })
        ));
        bad = p;
        bad.budget.maximum_free_dofs = 23;
        assert!(matches!(
            bad.assemble(cx),
            Err(TetElasticError::BudgetExceeded {
                what: "free_dofs",
                ..
            })
        ));
        bad = p;
        bad.fixed_dofs = &[1, 1];
        assert!(matches!(
            bad.assemble(cx),
            Err(TetElasticError::DuplicateFixedDof { dof: 1 })
        ));
        bad.fixed_dofs = &[24];
        assert!(matches!(
            bad.assemble(cx),
            Err(TetElasticError::FixedDofOutOfRange { .. })
        ));
        let fixed: Vec<_> = (0..24).collect();
        bad.fixed_dofs = &fixed;
        assert!(matches!(bad.assemble(cx), Err(TetElasticError::NoFreeDofs)));
        bad = p;
        bad.materials = TetMaterialField::PerElement(&[]);
        assert!(matches!(
            bad.assemble(cx),
            Err(TetElasticError::MaterialCountMismatch { .. })
        ));
        let repeated = [[0, 1, 2, 3, 4, 5, 6, 6]];
        bad = p;
        bad.hexahedra = &repeated;
        assert!(matches!(
            bad.assemble(cx),
            Err(TetElasticError::RepeatedTetNode { .. })
        ));
        let outside = [[0, 1, 2, 3, 4, 5, 6, 8]];
        bad.hexahedra = &outside;
        assert!(matches!(
            bad.assemble(cx),
            Err(TetElasticError::NodeIndexOutOfRange { node: 8, .. })
        ));
        let duplicate = [CELLS[0], CELLS[0]];
        bad.hexahedra = &duplicate;
        assert!(matches!(
            bad.assemble(cx),
            Err(TetElasticError::InvalidGeometryUpdate { .. })
        ));

        // All 2x2x2 Gauss Jacobians are positive, but the moved corner folds.
        let mut folded = CORNERS;
        folded[6] = [0.0; 3];
        let g = 1.0 / 3.0_f64.sqrt();
        for x in [-g, g] {
            for y in [-g, g] {
                for z in [-g, g] {
                    assert!(hex8_point(folded, [x, y, z], 1e-12, 0, 0).is_ok());
                }
            }
        }
        assert!(matches!(
            problem(&folded, &material, Formulation::BBar).assemble(cx),
            Err(TetElasticError::NonPositiveJacobian { .. })
        ));
    });
    with_cx(true, |cx| {
        assert!(matches!(
            problem(&CORNERS, &coupled_material(), Formulation::BBar).assemble(cx),
            Err(TetElasticError::Cancelled)
        ));
    });
}
