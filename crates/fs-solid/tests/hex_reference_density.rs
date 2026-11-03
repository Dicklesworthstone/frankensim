//! G0/G1/G3: graded physical inertia through both native quadratic hex families.
//! Exact monomial integrals are independent of the assembly quadrature.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solid::linear3::hex::{
    HEX20_REFERENCE_NODES, HEX27_REFERENCE_NODES, HexConnectivity, HexElasticProblem, HexFamily,
    HexQuadrature, hex20_shapes, hex27_shapes,
};
use fs_solid::linear3::{
    QuadraticThermalStrainField, TetAssemblyBudget, TetElasticAssembly, TetElasticError,
    TetElasticMaterial, TetMaterialField, TetThermalStrainState,
};

const FAMILIES: [HexFamily; 2] = [HexFamily::Serendipity20, HexFamily::Lagrange27];
const ELEMENTS20: [[usize; 20]; 1] = [[
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
]];
const ELEMENTS27: [[usize; 27]; 1] = [[
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26,
]];

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 61,
                    kernel_id: 61,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn local_nodes(family: HexFamily) -> &'static [[i8; 3]] {
    match family {
        HexFamily::Serendipity20 => &HEX20_REFERENCE_NODES,
        HexFamily::Lagrange27 => &HEX27_REFERENCE_NODES,
    }
}

fn problem<'a>(
    family: HexFamily,
    nodes: &'a [[f64; 3]],
    material: &'a TetElasticMaterial,
    fixed: &'a [usize],
    stiffness_points: usize,
) -> HexElasticProblem<'a> {
    HexElasticProblem {
        nodes_m: nodes,
        family,
        elements: match family {
            HexFamily::Serendipity20 => HexConnectivity::Hex20(&ELEMENTS20),
            HexFamily::Lagrange27 => HexConnectivity::Hex27(&ELEMENTS27),
        },
        materials: TetMaterialField::Uniform(material),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::standard(),
        quadrature: HexQuadrature {
            points_per_dir: stiffness_points,
        },
    }
}

fn material() -> TetElasticMaterial {
    TetElasticMaterial::try_new(2000.0, 12e9, 0.25, [7; 32]).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

type Polynomial = Vec<([usize; 3], f64)>;

fn multiply(a: &Polynomial, b: &Polynomial) -> Polynomial {
    a.iter()
        .flat_map(|(p, c)| {
            b.iter()
                .map(move |(q, d)| (core::array::from_fn(|i| p[i] + q[i]), c * d))
        })
        .collect()
}

fn integrate(poly: &Polynomial) -> f64 {
    poly.iter()
        .map(|(power, coefficient)| {
            coefficient
                * power
                    .iter()
                    .map(|p| {
                        if p % 2 == 0 {
                            2.0 / (p + 1) as f64
                        } else {
                            0.0
                        }
                    })
                    .product::<f64>()
        })
        .sum()
}

fn axis_polynomial(axis: usize, coefficients: &[f64]) -> Polynomial {
    coefficients
        .iter()
        .enumerate()
        .filter(|(_, c)| **c != 0.0)
        .map(|(degree, &c)| {
            let mut power = [0; 3];
            power[axis] = degree;
            (power, c)
        })
        .collect()
}

fn shape_polynomials(family: HexFamily) -> Vec<Polynomial> {
    local_nodes(family)
        .iter()
        .enumerate()
        .map(|(node, coordinate)| {
            if family == HexFamily::Lagrange27 {
                (0..3).fold(vec![([0; 3], 1.0)], |p, axis| {
                    let coefficients = match coordinate[axis] {
                        -1 => [0.0, -0.5, 0.5],
                        0 => [1.0, 0.0, -1.0],
                        1 => [0.0, 0.5, 0.5],
                        _ => unreachable!(),
                    };
                    multiply(&p, &axis_polynomial(axis, &coefficients))
                })
            } else if node < 8 {
                let product = (0..3).fold(vec![([0; 3], 0.125)], |p, axis| {
                    multiply(&p, &axis_polynomial(axis, &[1.0, coordinate[axis] as f64]))
                });
                let mut linear = vec![([0; 3], -2.0)];
                for axis in 0..3 {
                    linear.extend(axis_polynomial(axis, &[0.0, coordinate[axis] as f64]));
                }
                multiply(&product, &linear)
            } else {
                (0..3).fold(vec![([0; 3], 0.25)], |p, axis| {
                    let coefficients = if coordinate[axis] == 0 {
                        [1.0, 0.0, -1.0]
                    } else {
                        [1.0, coordinate[axis] as f64, 0.0]
                    };
                    multiply(&p, &axis_polynomial(axis, &coefficients))
                })
            }
        })
        .collect()
}

fn density(p: [f64; 3]) -> f64 {
    2000.0 + 500.0 * p[0] + 250.0 * p[1] + 125.0 * p[2] + 62.5 * p[0] * p[1] * p[2]
}

fn corner_density() -> [[f64; 8]; 1] {
    [core::array::from_fn(|i| {
        density(HEX27_REFERENCE_NODES[i].map(f64::from))
    })]
}

fn density_polynomial() -> Polynomial {
    vec![
        ([0, 0, 0], 2000.0),
        ([1, 0, 0], 500.0),
        ([0, 1, 0], 250.0),
        ([0, 0, 1], 125.0),
        ([1, 1, 1], 62.5),
    ]
}

fn physical_polynomials(curved: bool) -> [Polynomial; 3] {
    if curved {
        // These maps belong to BOTH serendipity20 and tensor Q2 spaces.
        [
            vec![([1, 0, 0], 1.0), ([2, 0, 0], 0.125)],
            vec![([0, 1, 0], 1.0), ([2, 1, 0], 0.25)],
            vec![([0, 0, 1], 1.0), ([2, 0, 1], 0.375)],
        ]
    } else {
        core::array::from_fn(|axis| axis_polynomial(axis, &[0.0, 1.0]))
    }
}

fn geometry(family: HexFamily, curved: bool) -> Vec<[f64; 3]> {
    local_nodes(family)
        .iter()
        .map(|p| {
            let [r, s, t] = p.map(f64::from);
            if curved {
                [
                    r + 0.125 * r * r,
                    s * (1.0 + 0.25 * r * r),
                    t * (1.0 + 0.375 * r * r),
                ]
            } else {
                [r, s, t]
            }
        })
        .collect()
}

fn measure(curved: bool) -> Polynomial {
    if curved {
        // J is triangular: detJ=(1+r/4)(1+r²/4)(1+3r²/8), degree FIVE.
        let jacobian = multiply(
            &multiply(
                &axis_polynomial(0, &[1.0, 0.25]),
                &axis_polynomial(0, &[1.0, 0.0, 0.25]),
            ),
            &axis_polynomial(0, &[1.0, 0.0, 0.375]),
        );
        multiply(&density_polynomial(), &jacobian)
    } else {
        density_polynomial()
    }
}

fn check_mass(family: HexFamily, nodes: &[[f64; 3]], a: &TetElasticAssembly, curved: bool) {
    let polynomials = shape_polynomials(family);
    let measure = measure(curved);
    let mass = a.mass.to_dense();
    let dofs = 3 * nodes.len();
    for (i, ni) in polynomials.iter().enumerate() {
        for (j, nj) in polynomials.iter().enumerate() {
            let exact = integrate(&multiply(&multiply(ni, nj), &measure));
            for c in 0..3 {
                for d in 0..3 {
                    close(
                        mass[(3 * i + c) * dofs + 3 * j + d],
                        if c == d { exact } else { 0.0 },
                    );
                }
            }
        }
    }
    close(a.total_mass_kg, integrate(&measure));
    let translation: Vec<_> = (0..dofs)
        .map(|i| if i % 3 == 0 { 1.0 } else { 0.0 })
        .collect();
    let mut momentum = vec![0.0; dofs];
    a.mass.spmv(&translation, &mut momentum);
    let physical = physical_polynomials(curved);
    for axis in 0..3 {
        let moment: f64 = nodes
            .iter()
            .enumerate()
            .map(|(i, p)| p[axis] * momentum[3 * i])
            .sum();
        close(moment, integrate(&multiply(&physical[axis], &measure)));
    }
    let rotation: Vec<_> = nodes.iter().flat_map(|p| [0.0, -p[2], p[1]]).collect();
    a.mass.spmv(&rotation, &mut momentum);
    let energy_twice: f64 = rotation.iter().zip(momentum).map(|(v, p)| v * p).sum();
    let mut radius_squared = multiply(&physical[1], &physical[1]);
    radius_squared.extend(multiply(&physical[2], &physical[2]));
    close(
        energy_twice,
        integrate(&multiply(&radius_squared, &measure)),
    );
    if family == HexFamily::Serendipity20 {
        assert!(
            a.mass.row(0).1.iter().sum::<f64>() < 0.0,
            "corner row sums must not be lumped"
        );
    }
}

#[test]
fn affine_q1_density_preserves_all_mass_entries_and_rigid_moments() {
    with_cx(|cx, _| {
        let material = material();
        for family in FAMILIES {
            let nodes = geometry(family, false);
            let p = problem(family, &nodes, &material, &[], 3);
            let a = p
                .assemble_with_reference_density(
                    &corner_density(),
                    HexQuadrature { points_per_dir: 3 },
                    cx,
                )
                .unwrap();
            check_mass(family, &nodes, &a, false);
            assert_eq!(
                a.stiffness.to_dense(),
                p.assemble(cx).unwrap().stiffness.to_dense()
            );
        }
    });
}

#[test]
fn curved_degree_ten_mass_uses_separate_exact_rule_without_extra_stiffness_work() {
    with_cx(|cx, _| {
        let material = material();
        for family in FAMILIES {
            let nodes = geometry(family, true);
            let p = problem(family, &nodes, &material, &[], 3);
            for points in [3, 4, 5] {
                assert!(matches!(
                    p.assemble_with_reference_density(
                        &corner_density(),
                        HexQuadrature {
                            points_per_dir: points
                        },
                        cx
                    ),
                    Err(TetElasticError::QuadratureInadmissible { .. })
                ));
            }
            let a = p
                .assemble_with_reference_density(
                    &corner_density(),
                    HexQuadrature { points_per_dir: 6 },
                    cx,
                )
                .unwrap();
            check_mass(family, &nodes, &a, true);
            assert_eq!(
                a.stiffness.to_dense(),
                p.assemble(cx).unwrap().stiffness.to_dense()
            );
            // The degree-ten edge mass really needs six points: five would
            // preserve total mass but miss this independently known entry.
            let under: f64 = HexQuadrature { points_per_dir: 5 }
                .points()
                .unwrap()
                .iter()
                .map(|(r, w)| {
                    let shape = match family {
                        HexFamily::Serendipity20 => hex20_shapes(*r).0[8],
                        HexFamily::Lagrange27 => hex27_shapes(*r).0[8],
                    };
                    let x = r[0];
                    w * shape
                        * shape
                        * density(*r)
                        * (1.0 + 0.25 * x)
                        * (1.0 + 0.25 * x * x)
                        * (1.0 + 0.375 * x * x)
                })
                .sum();
            assert!((under - a.mass.to_dense()[24 * (3 * nodes.len()) + 24]).abs() > 1e-5);
        }
    });
}

#[test]
fn density_placement_changes_native_supported_eigenvalue_without_changing_total_mass() {
    with_cx(|cx, _| {
        let material = material();
        let densities = [core::array::from_fn(|i| {
            2000.0 * (1.0 + 0.5 * f64::from(HEX27_REFERENCE_NODES[i][2]))
        })];
        for family in FAMILIES {
            let nodes = geometry(family, false);
            // One moving displacement at edge 8 gives the exact native
            // generalized eigenvalue K_ii/M_ii; no eigensolver approximation.
            let fixed: Vec<_> = (0..3 * nodes.len()).filter(|i| *i != 24).collect();
            let p = problem(family, &nodes, &material, &fixed, 3);
            let uniform = p.assemble(cx).unwrap();
            let graded = p
                .assemble_with_reference_density(
                    &densities,
                    HexQuadrature { points_per_dir: 3 },
                    cx,
                )
                .unwrap();
            assert_eq!(graded.free_dofs, [24]);
            assert_eq!(graded.stiffness.to_dense(), uniform.stiffness.to_dense());
            close(graded.total_mass_kg, uniform.total_mass_kg);
            let lambda = |a: &TetElasticAssembly| a.stiffness.to_dense()[0] / a.mass.to_dense()[0];
            let ratio = match family {
                HexFamily::Serendipity20 => 4.0 / 3.0,
                HexFamily::Lagrange27 => 8.0 / 5.0,
            };
            close(lambda(&graded) / lambda(&uniform), ratio);
        }
    });
}

#[test]
fn constant_profile_preserves_bits_but_a_different_mass_rule_is_honored() {
    with_cx(|cx, _| {
        let material = material();
        let density = [[material.density_kg_m3; 8]];
        for family in FAMILIES {
            for curved in [false, true] {
                let nodes = geometry(family, curved);
                let p = problem(family, &nodes, &material, &[], 3);
                let original = p.assemble(cx).unwrap();
                let unchanged = p
                    .assemble_with_reference_density(
                        &density,
                        HexQuadrature { points_per_dir: 3 },
                        cx,
                    )
                    .unwrap();
                assert_eq!(original.mass.to_dense(), unchanged.mass.to_dense());
                assert_eq!(
                    original.stiffness.to_dense(),
                    unchanged.stiffness.to_dense()
                );
                assert_eq!(
                    original.total_mass_kg.to_bits(),
                    unchanged.total_mass_kg.to_bits()
                );
                let refined = p
                    .assemble_with_reference_density(
                        &density,
                        HexQuadrature { points_per_dir: 6 },
                        cx,
                    )
                    .unwrap();
                let both_refined = problem(family, &nodes, &material, &[], 6)
                    .assemble(cx)
                    .unwrap();
                assert_eq!(refined.mass.to_dense(), both_refined.mass.to_dense());
                assert_eq!(refined.stiffness.to_dense(), original.stiffness.to_dense());
                if curved {
                    assert_ne!(refined.mass.to_dense(), original.mass.to_dense());
                }
            }
        }
    });
}

#[test]
fn frozen_graded_mass_uses_original_inventory_with_local_hot_stiffness() {
    with_cx(|cx, _| {
        let material = material();
        let hot = TetElasticMaterial::try_new(2000.0, 6e9, 0.25, [8; 32]).unwrap();
        let state = TetThermalStrainState::try_new(
            600.0,
            300.0,
            [0.001; 6],
            [8; 32],
            ContentHash([9; 32]),
            0.01,
        )
        .unwrap();
        let states = vec![state; 27];
        let materials = vec![hot; 27];
        let field = QuadraticThermalStrainField::ThermomechanicalAtQuadrature {
            states: &states,
            materials: &materials,
            points_per_element: 27,
        };
        for family in FAMILIES {
            let nodes = geometry(family, true);
            let p = problem(family, &nodes, &material, &[], 3);
            let mass_rule = HexQuadrature { points_per_dir: 6 };
            let cold = p
                .assemble_with_reference_density(&corner_density(), mass_rule, cx)
                .unwrap();
            let warmed = p
                .assemble_thermoelastic_with_reference_density(
                    field,
                    &corner_density(),
                    mass_rule,
                    cx,
                )
                .unwrap();
            assert_eq!(cold.mass.to_dense(), warmed.mass.to_dense());
            assert_eq!(cold.total_mass_kg.to_bits(), warmed.total_mass_kg.to_bits());
            for (before, after) in cold
                .stiffness
                .to_dense()
                .iter()
                .zip(warmed.stiffness.to_dense())
            {
                close(0.5 * before, after);
            }
        }
    });
}

#[test]
fn density_quadrature_budget_and_cancellation_refuse_without_a_partial_pencil() {
    for family in FAMILIES {
        with_cx(|cx, gate| {
            let material = material();
            let nodes = geometry(family, false);
            let mut p = problem(family, &nodes, &material, &[], 3);
            let q = HexQuadrature { points_per_dir: 3 };
            assert!(matches!(
                p.assemble_with_reference_density(&[], q, cx),
                Err(TetElasticError::MaterialCountMismatch { .. })
            ));
            for value in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
                assert!(matches!(
                    p.assemble_with_reference_density(&[[value; 8]], q, cx),
                    Err(TetElasticError::InvalidMaterial { .. })
                ));
            }
            assert!(matches!(
                p.assemble_with_reference_density(
                    &corner_density(),
                    HexQuadrature { points_per_dir: 2 },
                    cx
                ),
                Err(TetElasticError::UnderIntegrated { .. })
            ));
            assert!(matches!(
                p.assemble_with_reference_density(
                    &corner_density(),
                    HexQuadrature { points_per_dir: 13 },
                    cx
                ),
                Err(TetElasticError::QuadratureInadmissible { .. })
            ));
            p.budget.maximum_nodes = nodes.len() - 1;
            assert!(matches!(
                p.assemble_with_reference_density(&corner_density(), q, cx),
                Err(TetElasticError::BudgetExceeded { what: "nodes", .. })
            ));
            gate.request();
            assert!(matches!(
                p.assemble_with_reference_density(&corner_density(), q, cx),
                Err(TetElasticError::Cancelled)
            ));
        });
    }
}
