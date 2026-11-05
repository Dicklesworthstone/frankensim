//! Focused G0/G1/G3 shell spatial-operator battery (frankensim-b8bxd.9.1).
//!
//! Single-triangle checks know the six rigid modes without an eigen-solver
//! oracle. Affine patches on refined and rotated meshes independently exercise
//! plane-stress elasticity. Bending remains a surrogate with no continuum claim.

use fs_solid::{
    AssemblyBudget, DampingModel, OperatorDiagnostics, ShellError, ShellIdentity, ShellMaterial,
    ShellNode, ShellPlate, ShellSupport,
};

fn plate() -> ShellPlate {
    ShellPlate {
        nodes: vec![
            ShellNode {
                position_m: [0.0, 0.0, 0.0],
            },
            ShellNode {
                position_m: [1.0, 0.0, 0.0],
            },
            ShellNode {
                position_m: [0.0, 1.0, 0.0],
            },
        ],
        triangles: vec![[0, 1, 2]],
        thickness_m: 0.02,
        material: ShellMaterial {
            youngs_modulus_pa: 70.0e9,
            poisson_ratio: 0.33,
            density_kg_m3: 2_700.0,
        },
        identity: ShellIdentity {
            model_id: "euler-disc-base-v1".into(),
            source_id: "unit-flat-triangle".into(),
            state_id: "initial".into(),
        },
        support: None,
        damping: DampingModel::None,
        budget: AssemblyBudget::default(),
    }
}

fn rigid_displacement(plate: &ShellPlate, translation: [f64; 3], rotation: [f64; 3]) -> Vec<f64> {
    plate
        .nodes
        .iter()
        .flat_map(|node| {
            let [x, y, z] = node.position_m;
            let [rx, ry, rz] = rotation;
            [
                translation[0] + ry * z - rz * y,
                translation[1] + rz * x - rx * z,
                translation[2] + rx * y - ry * x,
                rx,
                ry,
                rz,
            ]
        })
        .collect()
}

fn rotated_plate() -> ShellPlate {
    let mut model = plate();
    let axis = [
        1.0 / 14.0_f64.sqrt(),
        -2.0 / 14.0_f64.sqrt(),
        3.0 / 14.0_f64.sqrt(),
    ];
    let angle = 0.73_f64;
    let (sine, cosine) = angle.sin_cos();
    for node in &mut model.nodes {
        let point = node.position_m;
        let axis_cross_point = [
            axis[1] * point[2] - axis[2] * point[1],
            axis[2] * point[0] - axis[0] * point[2],
            axis[0] * point[1] - axis[1] * point[0],
        ];
        let axis_dot_point = axis[0] * point[0] + axis[1] * point[1] + axis[2] * point[2];
        node.position_m = [
            cosine * point[0]
                + sine * axis_cross_point[0]
                + (1.0 - cosine) * axis_dot_point * axis[0]
                + 0.4,
            cosine * point[1]
                + sine * axis_cross_point[1]
                + (1.0 - cosine) * axis_dot_point * axis[1]
                - 0.7,
            cosine * point[2]
                + sine * axis_cross_point[2]
                + (1.0 - cosine) * axis_dot_point * axis[2]
                + 1.2,
        ];
    }
    model.identity.source_id = "arbitrarily-rotated-unit-flat-triangle".into();
    model
}

#[test]
fn shell_001_free_plate_has_typed_pd_inertia_six_rigid_modes_and_energy_consistency() {
    let model = plate();
    let assembly = model.assemble().expect("valid flat plate assembles");
    assert_eq!(assembly.full_mass.dimension(), 18);
    assert_eq!(assembly.full_stiffness.dimension(), 18);
    assert!(assembly.full_damping.is_none());
    assert_eq!(assembly.free_dofs, (0..18).collect::<Vec<_>>());
    let OperatorDiagnostics::Computed {
        raw_mass_min_eigenvalue,
        raw_stiffness_nullity,
        raw_stiffness_eigenvalue_spread,
        symmetry_residual,
    } = assembly.diagnostics
    else {
        panic!("fixture is within diagnostic budget");
    };
    assert!(
        raw_mass_min_eigenvalue > 0.0,
        "lumped mass array must be PD"
    );
    assert_eq!(
        raw_stiffness_nullity, 6,
        "single free panel has six rigid modes"
    );
    assert!(raw_stiffness_eigenvalue_spread.is_finite() && raw_stiffness_eigenvalue_spread > 1.0);
    assert!(symmetry_residual <= f64::EPSILON);

    // Formula-disjoint rigid-mode check: translations plus omega cross x,
    // rather than inferring the count from the diagnostic eigensweep.
    for (translation, rotation) in [
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
        ([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]),
        ([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
        ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        ([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        ([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
    ] {
        let vector = rigid_displacement(&model, translation, rotation);
        let energy = assembly.full_stiffness.quadratic_energy(&vector);
        assert!(
            energy.abs() < 1e-6,
            "rigid translation {translation:?}, rotation {rotation:?} has energy {energy:e}"
        );
    }

    let probe: Vec<f64> = (0..18).map(|index| (index as f64 - 7.0) * 0.03).collect();
    let work: f64 = probe
        .iter()
        .zip(assembly.full_stiffness.apply(&probe))
        .map(|(x, kx)| x * kx)
        .sum();
    assert!((2.0 * assembly.full_stiffness.quadratic_energy(&probe) - work).abs() < 1e-8);
    assert!(work >= -1e-8, "outer-product assembly is PSD");
}

#[test]
fn shell_001a_rotated_plate_preserves_all_global_rigid_modes() {
    let model = rotated_plate();
    let assembly = model.assemble().expect("rotated flat patch assembles");
    for (translation, rotation) in [
        ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        ([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        ([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
        ([0.7, -0.4, 0.2], [0.0, 0.0, 0.0]),
    ] {
        let energy = assembly
            .full_stiffness
            .quadratic_energy(&rigid_displacement(&model, translation, rotation));
        assert!(
            energy.abs() < 1.0e-4,
            "rotated rigid displacement has energy {energy:e}"
        );
    }
}

#[test]
fn shell_002_rayleigh_damping_and_three_point_pins_reduce_deterministically() {
    let mut model = plate();
    model.damping = DampingModel::Rayleigh {
        mass_proportional_per_s: 0.25,
        stiffness_proportional_s: 0.002,
    };
    model.support = Some(ShellSupport {
        node_indices: [2, 0, 1],
        normal: [0.0, 0.0, 1.0],
    });
    let assembly = model
        .assemble()
        .expect("supported Rayleigh plate assembles");
    assert_eq!(
        assembly.free_dofs,
        (3..18)
            .filter(|dof| ![6, 7, 8, 12, 13, 14].contains(dof))
            .collect::<Vec<_>>()
    );
    assert_eq!(assembly.mass.dimension(), 9);
    assert_eq!(assembly.stiffness.dimension(), 9);
    let damping = assembly
        .damping
        .as_ref()
        .expect("damping follows valid M and K");
    for (reduced_row, &full_row) in assembly.free_dofs.iter().enumerate() {
        for (reduced_column, &full_column) in assembly.free_dofs.iter().enumerate() {
            let reduced = damping.values()[reduced_row * damping.dimension() + reduced_column];
            let full = full_row * assembly.full_mass.dimension() + full_column;
            assert!(
                (reduced
                    - (0.25 * assembly.full_mass.values()[full]
                        + 0.002 * assembly.full_stiffness.values()[full]))
                    .abs()
                    < 1e-9
            );
        }
    }
    // Repeated assembly preserves both the declared element ordering and the
    // canonical node-major reduced DOF ordering.
    assert_eq!(assembly, model.assemble().expect("deterministic repeat"));
}

#[test]
fn shell_003_refuses_malformed_budgeted_and_out_of_applicability_requests() {
    let mut malformed = plate();
    malformed.material.poisson_ratio = 0.5;
    assert!(matches!(
        malformed.assemble(),
        Err(ShellError::InvalidInput { .. })
    ));

    let mut bad_support = plate();
    bad_support.support = Some(ShellSupport {
        node_indices: [0, 0, 2],
        normal: [0.0, 0.0, 1.0],
    });
    assert!(matches!(
        bad_support.assemble(),
        Err(ShellError::UnsupportedBoundary { .. })
    ));

    let mut opposite_support = plate();
    opposite_support.support = Some(ShellSupport {
        node_indices: [0, 1, 2],
        normal: [0.0, 0.0, -1.0],
    });
    assert!(matches!(
        opposite_support.assemble(),
        Err(ShellError::UnsupportedBoundary { .. })
    ));

    let mut unused_node = plate();
    unused_node.nodes.push(ShellNode {
        position_m: [0.25, 0.25, 0.0],
    });
    assert!(matches!(
        unused_node.assemble(),
        Err(ShellError::InvalidInput { .. })
    ));

    let mut duplicate = plate();
    duplicate.triangles.push([2, 1, 0]);
    assert!(matches!(
        duplicate.assemble(),
        Err(ShellError::InvalidInput { .. })
    ));

    let mut coincident = plate();
    coincident.nodes.push(ShellNode {
        position_m: [0.0, 0.0, 0.0],
    });
    coincident.triangles = vec![[0, 1, 3], [1, 2, 3], [2, 0, 3]];
    assert!(matches!(
        coincident.assemble(),
        Err(ShellError::InvalidInput { .. })
    ));

    let mut same_oriented_edge = plate();
    same_oriented_edge.nodes.push(ShellNode {
        position_m: [1.0, 1.0, 0.0],
    });
    same_oriented_edge.triangles.push([0, 1, 3]);
    assert!(matches!(
        same_oriented_edge.assemble(),
        Err(ShellError::UnsupportedGeometry { .. })
    ));

    let mut non_manifold = plate();
    non_manifold.nodes.extend([
        ShellNode {
            position_m: [1.0, 1.0, 0.0],
        },
        ShellNode {
            position_m: [0.5, -1.0, 0.0],
        },
    ]);
    non_manifold.triangles.extend([[1, 0, 3], [1, 0, 4]]);
    assert!(matches!(
        non_manifold.assemble(),
        Err(ShellError::UnsupportedGeometry { .. })
    ));

    let mut disconnected = plate();
    disconnected.nodes.extend([
        ShellNode {
            position_m: [2.0, 0.0, 0.0],
        },
        ShellNode {
            position_m: [3.0, 0.0, 0.0],
        },
        ShellNode {
            position_m: [2.0, 1.0, 0.0],
        },
    ]);
    disconnected.triangles.push([3, 4, 5]);
    assert!(matches!(
        disconnected.assemble(),
        Err(ShellError::UnsupportedGeometry { .. })
    ));

    let mut collinear_support = plate();
    collinear_support.nodes.push(ShellNode {
        position_m: [0.5, 0.0, 0.0],
    });
    collinear_support.triangles = vec![[0, 3, 2], [3, 1, 2]];
    collinear_support.support = Some(ShellSupport {
        node_indices: [0, 3, 1],
        normal: [0.0, 0.0, 1.0],
    });
    assert!(matches!(
        collinear_support.assemble(),
        Err(ShellError::UnsupportedBoundary { .. })
    ));

    let mut non_finite_derived = plate();
    non_finite_derived.thickness_m = 1.0e200;
    assert!(matches!(
        non_finite_derived.assemble(),
        Err(ShellError::InvalidInput { .. })
    ));

    let mut over_budget = plate();
    over_budget.budget.max_matrix_entries = 17 * 18;
    assert!(matches!(
        over_budget.assemble(),
        Err(ShellError::BudgetExceeded { .. })
    ));

    let mut validation_over_budget = plate();
    validation_over_budget.budget.max_work_units = 1;
    assert!(matches!(
        validation_over_budget.assemble(),
        Err(ShellError::BudgetExceeded { .. })
    ));

    // Two non-coplanar panels deliberately request a curved/folded shell;
    // this flat spatial slice must refuse rather than silently flatten it.
    let mut folded = plate();
    folded.nodes.push(ShellNode {
        position_m: [0.0, 0.0, 1.0],
    });
    folded.triangles.push([0, 3, 1]);
    assert!(matches!(
        folded.assemble(),
        Err(ShellError::UnsupportedGeometry { .. })
    ));
}

#[test]
fn shell_004_conditioning_can_be_explicitly_bounded_without_refusing_assembly() {
    let mut model = plate();
    model.budget.max_conditioning_dofs = 17;
    let assembly = model
        .assemble()
        .expect("conditioning is optional within its explicit budget");
    assert!(matches!(
        assembly.diagnostics,
        OperatorDiagnostics::NotComputed { .. }
    ));
}

#[test]
fn shell_g1_membrane_affine_patch_has_exact_plane_stress_energy_under_refinement_and_rotation() {
    // Formula-disjoint oracle: integrate the continuum plane-stress density
    // over a parallelogram, independent of triangle shapes or stiffness entries.
    let frames = [
        ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
        (
            [1.0 / 3.0, 2.0 / 3.0, -2.0 / 3.0],
            [2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0],
        ),
    ];
    for n in [1, 2, 4] {
        for (tangent, bitangent) in frames {
            for poisson in [-0.7, 0.0, 0.33, 0.499] {
                let mut model = plate();
                model.material.poisson_ratio = poisson;
                model.budget.max_conditioning_dofs = 0;
                let coordinates: Vec<[f64; 2]> = (0..=n)
                    .flat_map(|j| {
                        (0..=n).map(move |i| {
                            [
                                (i as f64 + 0.3 * j as f64) / n as f64,
                                0.8 * j as f64 / n as f64,
                            ]
                        })
                    })
                    .collect();
                model.nodes = coordinates
                    .iter()
                    .map(|&[x, y]| ShellNode {
                        position_m: core::array::from_fn(|c| {
                            [0.4, -0.7, 1.2][c] + tangent[c] * x + bitangent[c] * y
                        }),
                    })
                    .collect();
                model.triangles = (0..n)
                    .flat_map(|j| {
                        (0..n).flat_map(move |i| {
                            let a = j * (n + 1) + i;
                            [[a, a + 1, a + n + 2], [a, a + n + 2, a + n + 1]]
                        })
                    })
                    .collect();
                let stiffness = model.assemble().unwrap().full_stiffness;
                for [ex, ey, exy] in [
                    [0.002, -poisson * 0.002, 0.0],
                    [0.002, 0.002, 0.0],
                    [0.0, 0.0, 0.003],
                    [0.002, -0.001, 0.003],
                ] {
                    let u: Vec<_> = coordinates
                        .iter()
                        .flat_map(|&[x, y]| {
                            let ut = ex * x + exy * y;
                            let ub = exy * x + ey * y;
                            [
                                ut * tangent[0] + ub * bitangent[0],
                                ut * tangent[1] + ub * bitangent[1],
                                ut * tangent[2] + ub * bitangent[2],
                                0.0,
                                0.0,
                                0.0,
                            ]
                        })
                        .collect();
                    let e = model.material.youngs_modulus_pa;
                    let expected = 0.8
                        * model.thickness_m
                        * e
                        * (0.5 * (ex * ex + ey * ey + 2.0 * poisson * ex * ey)
                            / (1.0 - poisson * poisson)
                            + exy * exy / (1.0 + poisson));
                    let energy = stiffness.quadratic_energy(&u);
                    assert!(
                        (energy / expected - 1.0).abs() < 2e-12,
                        "n={n}, nu={poisson}, strain={:?}: {energy} != {expected}",
                        [ex, ey, exy]
                    );
                    let forces = stiffness.apply(&u);
                    let force_scale = e * model.thickness_m * 0.003;
                    for j in 1..n {
                        for i in 1..n {
                            let node = j * (n + 1) + i;
                            assert!(
                                forces[6 * node..6 * node + 6]
                                    .iter()
                                    .all(|f| f.abs() < 2e-12 * force_scale),
                                "constant-stress patch has an unbalanced interior node"
                            );
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn shell_g1_membrane_uniaxial_traction_recovers_poisson_contraction() {
    for poisson in [-0.7, 0.0, 0.33, 0.499] {
        let mut model = plate();
        model.material.poisson_ratio = poisson;
        let strain = 0.002;
        let u: Vec<_> = model
            .nodes
            .iter()
            .flat_map(|node| {
                [
                    strain * node.position_m[0],
                    -poisson * strain * node.position_m[1],
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ]
            })
            .collect();
        let force = model.assemble().unwrap().full_stiffness.apply(&u);
        let axial = 0.5 * model.material.youngs_modulus_pa * model.thickness_m * strain;
        for (index, actual) in force.into_iter().enumerate() {
            let expected = match index {
                0 => -axial,
                6 => axial,
                _ => 0.0,
            };
            assert!(
                (actual - expected).abs() < 1e-12 * axial,
                "nu={poisson}, dof={index}: {actual} != {expected}"
            );
        }
    }
}
