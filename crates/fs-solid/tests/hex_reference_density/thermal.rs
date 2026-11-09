//! Original Q1 inertia through solved thermal/dead-load/support prestress.

use super::*;
use fs_solid::linear::Formulation;
use fs_solid::linear3::hex8::Hex8ElasticProblem;
use fs_solid::linear3::{TetStaticSolveConfig, TetThermalStrainField};

const FRAME: ContentHash = ContentHash([62; 32]);
const CELLS8: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];

fn thermal_state(law: &TetElasticMaterial, alpha: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        300.0 + 1000.0 * alpha,
        300.0,
        [alpha, alpha, alpha, 0.0, 0.0, 0.0],
        law.material_state_identity,
        ContentHash([63; 32]),
        0.02,
    )
    .unwrap()
}

fn fixed_except_first(nodes: usize) -> Vec<usize> {
    (1..3 * nodes).collect()
}

fn derivative(p: &Polynomial, axis: usize) -> Polynomial {
    p.iter()
        .filter(|(power, _)| power[axis] != 0)
        .map(|(power, coefficient)| {
            let mut next = *power;
            next[axis] -= 1;
            (next, coefficient * power[axis] as f64)
        })
        .collect()
}

fn exact_compressive_stiffness(shape: &Polynomial, alpha: f64) -> f64 {
    // E=12e9 and nu=1/4 give sigma=-24e9*alpha I at zero displacement.
    -24e9
        * alpha
        * (0..3)
            .map(|axis| {
                let gradient = derivative(shape, axis);
                integrate(&multiply(&gradient, &gradient))
            })
            .sum::<f64>()
}

#[test]
fn graded_thermal_pencils_preserve_exact_mass_and_signed_prestress_moments() {
    with_cx(|cx, _| {
        let law = material();
        let alpha = 0.001;
        let state = thermal_state(&law, alpha);
        let field = TetThermalStrainField::Uniform(&state);
        let densities = corner_density();
        let rule = HexQuadrature { points_per_dir: 3 };
        for family in FAMILIES {
            let nodes = geometry(family, false);
            let fixed = fixed_except_first(nodes.len());
            let p = problem(family, &nodes, &law, &fixed, 3);
            // Dead load cancels the single free thermal load exactly. The
            // native equilibrium owner recovers uniform compressive stress.
            let forces = [(
                0,
                -p.assemble_thermal_load(field, cx).unwrap().reduced_force_n[0],
            )];
            let graded = p
                .solve_thermal_tangent_with_reference_density(
                    field,
                    &densities,
                    rule,
                    &forces,
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap();
            let scalar = p
                .solve_thermal_tangent_with_nodal_forces(
                    field,
                    &forces,
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap();
            let shape = &shape_polynomials(family)[0];
            let expected_mass = integrate(&multiply(&multiply(shape, shape), &measure(false)));
            close(graded.pencil().mass.get(0, 0), expected_mass);
            close(graded.pencil().total_mass_kg, scalar.pencil().total_mass_kg);
            assert!(
                (graded.pencil().mass.get(0, 0) / scalar.pencil().mass.get(0, 0) - 1.0).abs()
                    > 0.01
            );
            assert_eq!(graded.response(), scalar.response());
            assert_eq!(
                graded.pencil().stiffness.to_dense(),
                scalar.pencil().stiffness.to_dense()
            );
            close(
                graded.geometric_stiffness().get(0, 0),
                exact_compressive_stiffness(shape, alpha),
            );
            let lambda = graded.pencil().stiffness.get(0, 0) / graded.pencil().mass.get(0, 0);
            let scalar_lambda =
                scalar.pencil().stiffness.get(0, 0) / scalar.pencil().mass.get(0, 0);
            close(
                lambda / scalar_lambda,
                scalar.pencil().mass.get(0, 0) / expected_mass,
            );
        }
    });
}

#[test]
fn curved_thermal_pencils_keep_degree_ten_mass_and_original_stress_rule() {
    with_cx(|cx, _| {
        let law = material();
        let state = thermal_state(&law, 0.001);
        let field = TetThermalStrainField::Uniform(&state);
        for family in FAMILIES {
            let nodes = geometry(family, true);
            let fixed = fixed_except_first(nodes.len());
            let p = problem(family, &nodes, &law, &fixed, 3);
            let forces = [(
                0,
                -p.assemble_thermal_load(field, cx).unwrap().reduced_force_n[0],
            )];
            let solve = |points| {
                p.solve_thermal_tangent_with_reference_density(
                    field,
                    &corner_density(),
                    HexQuadrature {
                        points_per_dir: points,
                    },
                    &forces,
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.02,
                    cx,
                )
            };
            for points in [3, 4, 5] {
                assert!(matches!(
                    solve(points),
                    Err(TetElasticError::QuadratureInadmissible { .. })
                ));
            }
            let a = solve(6).unwrap();
            let b = solve(7).unwrap();
            let shape = &shape_polynomials(family)[0];
            let exact = integrate(&multiply(&multiply(shape, shape), &measure(true)));
            close(a.pencil().mass.get(0, 0), exact);
            close(b.pencil().mass.get(0, 0), exact);
            close(a.pencil().total_mass_kg, integrate(&measure(true)));
            assert_eq!(
                a.pencil().stiffness.to_dense(),
                b.pencil().stiffness.to_dense()
            );
            assert_eq!(
                a.geometric_stiffness().to_dense(),
                b.geometric_stiffness().to_dense()
            );
            assert_eq!(a.response(), b.response());
            assert_ne!(a.identity(), b.identity());
        }
    });
}

#[test]
fn prescribed_support_prestress_and_hot_materials_cannot_refresh_hex_inventory() {
    with_cx(|cx, _| {
        let laws = [
            material(),
            TetElasticMaterial::try_new(1200.0, 12e9, 0.25, [7; 32]).unwrap(),
            TetElasticMaterial::try_new(1200.0, 18e9, 0.25, [8; 32]).unwrap(),
        ];
        let forces = [(0, 12_000.0)];
        let prescribed = [(3, 0.002)];
        for family in FAMILIES {
            let nodes = geometry(family, false);
            let fixed = fixed_except_first(nodes.len());
            let mut masses = Vec::new();
            let mut stiffnesses = Vec::new();
            for law in &laws {
                let state = thermal_state(law, 0.001);
                let field = TetThermalStrainField::Uniform(&state);
                let p = problem(family, &nodes, law, &fixed, 3);
                let a = p
                    .solve_thermal_tangent_with_reference_density_and_prescribed_dofs(
                        field,
                        &corner_density(),
                        HexQuadrature { points_per_dir: 3 },
                        &forces,
                        &prescribed,
                        0.01,
                        TetStaticSolveConfig::default(),
                        FRAME,
                        0.02,
                        cx,
                    )
                    .unwrap();
                let scalar = p
                    .solve_thermal_tangent_with_prescribed_dofs(
                        field,
                        &forces,
                        &prescribed,
                        0.01,
                        TetStaticSolveConfig::default(),
                        FRAME,
                        0.02,
                        cx,
                    )
                    .unwrap();
                assert_eq!(a.response(), scalar.response());
                assert_eq!(
                    a.pencil().stiffness.to_dense(),
                    scalar.pencil().stiffness.to_dense()
                );
                assert_eq!(a.response().displacement_m()[1][0], 0.002);
                assert!(a.response().iterations() > 0);
                masses.push(a.pencil().mass.to_dense());
                stiffnesses.push(a.pencil().stiffness.get(0, 0));
                close(a.pencil().total_mass_kg, 16_000.0);
            }
            assert_eq!(masses[0], masses[1]);
            assert_eq!(masses[1], masses[2]);
            assert_eq!(stiffnesses[0], stiffnesses[1]);
            assert!(stiffnesses[2] > 1.4 * stiffnesses[1]);
        }
    });
}

fn hex8_problem<'a>(
    nodes: &'a [[f64; 3]],
    law: &'a TetElasticMaterial,
    fixed: &'a [usize],
    form: Formulation,
) -> Hex8ElasticProblem<'a> {
    Hex8ElasticProblem {
        nodes_m: nodes,
        hexahedra: &CELLS8,
        materials: TetMaterialField::Uniform(law),
        formulation: form,
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::standard(),
        stiffness_points_per_dir: 2,
    }
}

fn hex8_corner_shape() -> Polynomial {
    (0..3).fold(vec![([0; 3], 0.125)], |p, axis| {
        multiply(&p, &axis_polynomial(axis, &[1.0, -1.0]))
    })
}

#[test]
fn hex8_compatible_and_bbar_prestress_preserve_graded_affine_and_warped_mass() {
    with_cx(|cx, _| {
        let law = material();
        let alpha = 0.001;
        let state = thermal_state(&law, alpha);
        let field = TetThermalStrainField::Uniform(&state);
        let fixed = fixed_except_first(8);
        let shape = hex8_corner_shape();
        for warped in [false, true] {
            let nodes: Vec<_> = HEX27_REFERENCE_NODES[..8]
                .iter()
                .map(|p| {
                    let [r, s, t] = p.map(f64::from);
                    if warped {
                        [r + 0.125 * s * t, s + 0.25 * r * t, t + 0.375 * r * s]
                    } else {
                        [r, s, t]
                    }
                })
                .collect();
            let jacobian = if warped {
                // Exact determinant of the trilinear map, independent of native geometry.
                vec![
                    ([0; 3], 1.0),
                    ([2, 0, 0], -0.25 * 0.375),
                    ([0, 2, 0], -0.125 * 0.375),
                    ([0, 0, 2], -0.125 * 0.25),
                    ([1, 1, 1], 2.0 * 0.125 * 0.25 * 0.375),
                ]
            } else {
                vec![([0; 3], 1.0)]
            };
            let measure = multiply(&density_polynomial(), &jacobian);
            let exact = integrate(&multiply(&multiply(&shape, &shape), &measure));
            for form in [Formulation::Standard, Formulation::BBar] {
                let p = hex8_problem(&nodes, &law, &fixed, form);
                let forces = [(
                    0,
                    -p.assemble_thermal_load(field, cx).unwrap().reduced_force_n[0],
                )];
                let a = p
                    .solve_thermal_tangent_with_reference_density(
                        field,
                        &corner_density(),
                        &forces,
                        TetStaticSolveConfig::default(),
                        FRAME,
                        0.02,
                        cx,
                    )
                    .unwrap();
                let scalar = p
                    .solve_thermal_tangent_with_nodal_forces(
                        field,
                        &forces,
                        TetStaticSolveConfig::default(),
                        FRAME,
                        0.02,
                        cx,
                    )
                    .unwrap();
                close(a.pencil().mass.get(0, 0), exact);
                close(a.pencil().total_mass_kg, integrate(&measure));
                assert_eq!(a.response(), scalar.response());
                assert_eq!(
                    a.pencil().stiffness.to_dense(),
                    scalar.pencil().stiffness.to_dense()
                );
                assert_eq!(a.response().formulation(), form);
                if !warped {
                    close(
                        a.geometric_stiffness().get(0, 0),
                        exact_compressive_stiffness(&shape, alpha),
                    );
                    close(a.pencil().total_mass_kg, scalar.pencil().total_mass_kg);
                }
                assert!(
                    (a.pencil().mass.get(0, 0) / scalar.pencil().mass.get(0, 0) - 1.0).abs() > 0.01
                );
            }
        }
    });
}

#[test]
fn hex8_nonzero_support_and_material_updates_retain_original_inventory() {
    with_cx(|cx, _| {
        let nodes: Vec<_> = HEX27_REFERENCE_NODES[..8]
            .iter()
            .map(|p| p.map(f64::from))
            .collect();
        let fixed = fixed_except_first(8);
        let laws = [
            material(),
            TetElasticMaterial::try_new(1200.0, 18e9, 0.25, [8; 32]).unwrap(),
        ];
        for form in [Formulation::Standard, Formulation::BBar] {
            let mut masses = Vec::new();
            let mut stiffnesses = Vec::new();
            for law in &laws {
                let state = thermal_state(law, 0.001);
                let field = TetThermalStrainField::Uniform(&state);
                let p = hex8_problem(&nodes, law, &fixed, form);
                let a = p
                    .solve_thermal_tangent_with_reference_density_and_prescribed_dofs(
                        field,
                        &corner_density(),
                        &[(0, 12_000.0)],
                        &[(3, 0.002)],
                        0.01,
                        TetStaticSolveConfig::default(),
                        FRAME,
                        0.02,
                        cx,
                    )
                    .unwrap();
                let scalar = p
                    .solve_thermal_tangent_with_prescribed_dofs(
                        field,
                        &[(0, 12_000.0)],
                        &[(3, 0.002)],
                        0.01,
                        TetStaticSolveConfig::default(),
                        FRAME,
                        0.02,
                        cx,
                    )
                    .unwrap();
                assert_eq!(a.response(), scalar.response());
                assert_eq!(
                    a.pencil().stiffness.to_dense(),
                    scalar.pencil().stiffness.to_dense()
                );
                assert_eq!(a.response().displacement_m()[1][0], 0.002);
                assert!(a.response().iterations() > 0);
                masses.push(a.pencil().mass.to_dense());
                stiffnesses.push(a.pencil().stiffness.get(0, 0));
                close(a.pencil().total_mass_kg, 16_000.0);
            }
            assert_eq!(masses[0], masses[1]);
            assert!(stiffnesses[1] > 1.4 * stiffnesses[0]);
        }
    });
}

#[test]
fn density_is_admitted_before_equilibrium_and_scalar_operator_compatibility_is_retained() {
    with_cx(|cx, gate| {
        let law = material();
        let state = thermal_state(&law, 0.0);
        let field = TetThermalStrainField::Uniform(&state);
        let invalid_config = TetStaticSolveConfig {
            maximum_iterations: 0,
            ..Default::default()
        };
        for family in FAMILIES {
            let nodes = geometry(family, false);
            let fixed = fixed_except_first(nodes.len());
            let p = problem(family, &nodes, &law, &fixed, 3);
            let solve = |density: &[[f64; 8]], config| {
                p.solve_thermal_tangent_with_reference_density(
                    field,
                    density,
                    HexQuadrature { points_per_dir: 3 },
                    &[],
                    config,
                    FRAME,
                    0.02,
                    cx,
                )
            };
            assert!(matches!(
                solve(&[], invalid_config),
                Err(TetElasticError::MaterialCountMismatch { .. })
            ));
            for rho in [0.0, -1.0, f64::NAN, f64::INFINITY] {
                assert!(matches!(
                    solve(&[[rho; 8]], invalid_config),
                    Err(TetElasticError::InvalidMaterial { .. })
                ));
            }
            let a = solve(&[[2000.0; 8]], TetStaticSolveConfig::default()).unwrap();
            let b = solve(&[[2000.0; 8]], TetStaticSolveConfig::default()).unwrap();
            let scalar = p
                .solve_thermal_tangent_with_nodal_forces(
                    field,
                    &[],
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap();
            assert_eq!(a.pencil().mass.to_dense(), scalar.pencil().mass.to_dense());
            assert_eq!(
                a.pencil().stiffness.to_dense(),
                scalar.pencil().stiffness.to_dense()
            );
            assert_eq!(a.identity(), b.identity());
            assert_ne!(a.identity(), scalar.identity());
        }
        let nodes: Vec<_> = HEX27_REFERENCE_NODES[..8]
            .iter()
            .map(|p| p.map(f64::from))
            .collect();
        let fixed = fixed_except_first(8);
        let p = hex8_problem(&nodes, &law, &fixed, Formulation::BBar);
        let solve = |density: &[[f64; 8]], config| {
            p.solve_thermal_tangent_with_reference_density(
                field,
                density,
                &[],
                config,
                FRAME,
                0.02,
                cx,
            )
        };
        assert!(matches!(
            solve(&[], invalid_config),
            Err(TetElasticError::MaterialCountMismatch { .. })
        ));
        for rho in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                solve(&[[rho; 8]], invalid_config),
                Err(TetElasticError::InvalidMaterial { .. })
            ));
        }
        let a = solve(&[[2000.0; 8]], TetStaticSolveConfig::default()).unwrap();
        let scalar = p
            .solve_thermal_tangent_with_nodal_forces(
                field,
                &[],
                TetStaticSolveConfig::default(),
                FRAME,
                0.02,
                cx,
            )
            .unwrap();
        assert_eq!(a.pencil().mass.to_dense(), scalar.pencil().mass.to_dense());
        assert_eq!(
            a.pencil().stiffness.to_dense(),
            scalar.pencil().stiffness.to_dense()
        );
        assert_ne!(a.identity(), scalar.identity());
        gate.request();
        assert!(matches!(
            solve(&corner_density(), TetStaticSolveConfig::default()),
            Err(TetElasticError::Cancelled)
        ));
        for family in FAMILIES {
            let nodes = geometry(family, false);
            let fixed = fixed_except_first(nodes.len());
            let p = problem(family, &nodes, &law, &fixed, 3);
            assert!(matches!(
                p.solve_thermal_tangent_with_reference_density(
                    field,
                    &corner_density(),
                    HexQuadrature { points_per_dir: 3 },
                    &[],
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.02,
                    cx,
                ),
                Err(TetElasticError::Cancelled)
            ));
        }
    });
}

#[test]
fn supported_cells_keep_density_identity_when_every_reduced_mass_entry_is_hidden() {
    with_cx(|cx, _| {
        let law = material();
        let state = thermal_state(&law, 0.0);
        let field = TetThermalStrainField::Uniform(&state);
        let profile_a = [corner_density()[0], [2000.0; 8]];
        let profile_b = [corner_density()[0], corner_density()[0]];
        for family in FAMILIES {
            let first = geometry(family, false);
            let nodes: Vec<_> = first
                .iter()
                .copied()
                .chain(first.iter().map(|p| [p[0] + 4.0, p[1], p[2]]))
                .collect();
            let fixed = fixed_except_first(nodes.len());
            let cells20 = [ELEMENTS20[0], ELEMENTS20[0].map(|i| i + 20)];
            let cells27 = [ELEMENTS27[0], ELEMENTS27[0].map(|i| i + 27)];
            let p = HexElasticProblem {
                nodes_m: &nodes,
                family,
                elements: match family {
                    HexFamily::Serendipity20 => HexConnectivity::Hex20(&cells20),
                    HexFamily::Lagrange27 => HexConnectivity::Hex27(&cells27),
                },
                materials: TetMaterialField::Uniform(&law),
                fixed_dofs: &fixed,
                budget: TetAssemblyBudget::standard(),
                quadrature: HexQuadrature { points_per_dir: 3 },
            };
            let solve = |density: &[[f64; 8]]| {
                p.solve_thermal_tangent_with_reference_density(
                    field,
                    density,
                    HexQuadrature { points_per_dir: 3 },
                    &[],
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap()
            };
            let a = solve(&profile_a);
            let b = solve(&profile_b);
            assert_eq!(a.response(), b.response());
            assert_eq!(a.pencil().mass.to_dense(), b.pencil().mass.to_dense());
            assert_eq!(
                a.pencil().stiffness.to_dense(),
                b.pencil().stiffness.to_dense()
            );
            close(a.pencil().total_mass_kg, b.pencil().total_mass_kg);
            assert_ne!(
                a.identity(),
                b.identity(),
                "hidden inventory must remain identified"
            );
        }
        let first: Vec<_> = HEX27_REFERENCE_NODES[..8]
            .iter()
            .map(|p| p.map(f64::from))
            .collect();
        let nodes: Vec<_> = first
            .iter()
            .copied()
            .chain(first.iter().map(|p| [p[0] + 4.0, p[1], p[2]]))
            .collect();
        let fixed = fixed_except_first(nodes.len());
        let cells = [CELLS8[0], CELLS8[0].map(|i| i + 8)];
        for form in [Formulation::Standard, Formulation::BBar] {
            let p = Hex8ElasticProblem {
                hexahedra: &cells,
                ..hex8_problem(&nodes, &law, &fixed, form)
            };
            let solve = |density: &[[f64; 8]]| {
                p.solve_thermal_tangent_with_reference_density(
                    field,
                    density,
                    &[],
                    TetStaticSolveConfig::default(),
                    FRAME,
                    0.02,
                    cx,
                )
                .unwrap()
            };
            let a = solve(&profile_a);
            let b = solve(&profile_b);
            assert_eq!(a.response(), b.response());
            assert_eq!(a.pencil().mass.to_dense(), b.pencil().mass.to_dense());
            assert_eq!(
                a.pencil().stiffness.to_dense(),
                b.pencil().stiffness.to_dense()
            );
            close(a.pencil().total_mass_kg, b.pencil().total_mass_kg);
            assert_ne!(
                a.identity(),
                b.identity(),
                "hidden inventory must remain identified"
            );
        }
    });
}
