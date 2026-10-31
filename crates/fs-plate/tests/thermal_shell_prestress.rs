//! Planar initial-stress oracles on the native MITC pencil; synthetic sections.
use fs_material::sheet::SheetMatrix;
use fs_plate::PlateError;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::thermal::prestress::MembranePrestressLimits;
use fs_plate::shell::mitc::thermal::{FacetThermalStrain, ShellThermalStrain};
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellSection,
    ShellConstraint, ShellMass, assemble_mitc3p,
};
use fs_sparse::Csr;

const A: f64 = 1.0e6;
const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn limits() -> MembranePrestressLimits {
    MembranePrestressLimits {
        maximum_displacement_gradient: 0.01,
        maximum_rotation_rad: 0.01,
    }
}
fn section() -> ShearShellSection {
    ShearShellSection::new(
        SheetMatrix::orthotropic(A, A, 0.0, 0.5 * A).unwrap(),
        SheetMatrix::orthotropic(1.0, 1.0, 0.0, 0.5).unwrap(),
        [100.0, 0.0, 100.0],
        1.0,
        2.0,
        0.001,
        0.0,
    )
    .unwrap()
}
fn mesh() -> ShellMesh {
    ShellMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 0.0],
        ],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
    )
    .unwrap()
}
// Rows are support axes in world coordinates; the same frame maps the panel.
fn rotate(frame: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|k| (0..3).map(|j| frame[j][k] * v[j]).sum())
}
fn model(
    mesh: &ShellMesh,
    frame: [[f64; 3]; 3],
    normal_only: bool,
    center_only: bool,
) -> MitcShellModel {
    let directors = NodalDirectors::new(mesh, vec![frame[2]; mesh.nodes.len()], 0.9).unwrap();
    let constraints = if normal_only {
        (0..mesh.nodes.len())
            .map(|node| ShellConstraint {
                node,
                mask: [true, true, center_only && node != 4, true, true, true],
                frame: ConstraintFrame::Local(frame),
            })
            .collect()
    } else {
        vec![ShellConstraint {
            node: 0,
            mask: [true; 6],
            frame: ConstraintFrame::Local(frame),
        }]
    };
    assemble_mitc3p(
        mesh,
        &directors,
        &vec![section(); mesh.tris.len()],
        &constraints,
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: 5,
            max_triangles: 4,
        },
        1e-3,
        || false,
    )
    .unwrap()
}
fn uniform(free: [f64; 3], count: usize) -> Vec<FacetThermalStrain> {
    vec![ShellThermalStrain::try_new(free, [0.0; 3]).unwrap().into(); count]
}
fn entry(matrix: &Csr, row: usize, column: usize) -> f64 {
    let (columns, values) = matrix.row(row);
    columns
        .iter()
        .zip(values)
        .find_map(|(&c, &v)| (c == column).then_some(v))
        .unwrap_or(0.0)
}
fn action(matrix: &Csr, q: &[f64]) -> Vec<f64> {
    (0..q.len())
        .map(|row| {
            let (columns, values) = matrix.row(row);
            columns.iter().zip(values).map(|(&c, &v)| v * q[c]).sum()
        })
        .collect()
}
fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_constrained_heating_softens_the_native_mode_and_preserves_negative_stability_eigenvalues() {
    let mesh = mesh();
    let model = model(&mesh, IDENTITY, true, true);
    assert_eq!(model.pencil().free, 1);
    let k0 = entry(&model.pencil().k, 0, 0);
    let mass = entry(&model.pencil().m, 0, 0);
    let critical_free_strain = k0 / (4.0 * A);
    assert!(critical_free_strain > 0.0 && critical_free_strain < 0.001);
    for ratio in [-0.5, 0.0, 0.5, 2.0] {
        let free = ratio * critical_free_strain;
        let prepared = model
            .prepare_thermal(&uniform([free, free, 0.0], 4), || false)
            .unwrap();
        let tangent = prepared
            .solve_membrane_prestress(&[0.0; 30], limits(), || false)
            .unwrap();
        // Four quarter-area triangles have |grad(N_center)|=2. The exact
        // scalar eigenvalue is (K0 + 4*N)/M, with N=-A*integrated free strain.
        near(
            entry(tangent.geometric_stiffness(), 0, 0),
            -4.0 * A * free,
            3e-14,
            1e-12,
        );
        let lambda = entry(&tangent.pencil().k, 0, 0) / mass;
        near(lambda, (1.0 - ratio) * k0 / mass, 5e-14, 2e-11);
        if ratio == 2.0 {
            assert!(lambda < 0.0);
        } else if ratio == -0.5 {
            assert!(lambda > k0 / mass);
        }
        assert_eq!(tangent.pencil().m.row(0), model.pencil().m.row(0));
        assert_eq!(tangent.pencil().dof_map, model.pencil().dof_map);
        assert_eq!(tangent.elastic_stiffness().row(0), model.pencil().k.row(0));
        for facet in &tangent.reference_state().facets {
            near(facet.membrane_force_n_m[0], -A * free, 3e-14, 1e-12);
            near(facet.membrane_force_n_m[1], -A * free, 3e-14, 1e-12);
        }
        assert!(
            tangent
                .reference_state()
                .displacement_6n
                .iter()
                .all(|v| *v == 0.0)
        );
        assert_eq!(
            tangent.expand(&[0.003]).unwrap(),
            model.expand(&[0.003]).unwrap()
        );
    }
    near(entry(&model.pencil().k, 0, 0), k0, 0.0, 0.0);
}

#[test]
fn g1_free_expansion_has_no_membrane_geometric_stiffness_and_keeps_reference_inertia() {
    let mesh = mesh();
    let model = model(&mesh, IDENTITY, false, false);
    for free in [0.0001, -0.0001] {
        let prepared = model
            .prepare_thermal(&uniform([free, free, 0.0], 4), || false)
            .unwrap();
        let tangent = prepared
            .solve_membrane_prestress(&[0.0; 30], limits(), || false)
            .unwrap();
        for (node, p) in mesh.nodes.iter().enumerate() {
            near(
                tangent.reference_state().displacement_6n[6 * node],
                free * p[0],
                2e-8,
                1e-12,
            );
            near(
                tangent.reference_state().displacement_6n[6 * node + 1],
                free * p[1],
                2e-8,
                1e-12,
            );
        }
        for row in 0..model.pencil().free {
            for &value in tangent.geometric_stiffness().row(row).1 {
                near(value, 0.0, 0.0, 1e-6);
            }
            assert_eq!(tangent.pencil().m.row(row), model.pencil().m.row(row));
        }
    }
}

fn membrane_energy(mesh: &ShellMesh, normal: [f64; 3], full: &[f64], stress: [f64; 3]) -> f64 {
    mesh.tris
        .iter()
        .enumerate()
        .map(|(element, tri)| {
            let g = mesh.facet(element).unwrap();
            let w = tri.map(|node| (0..3).map(|j| normal[j] * full[6 * node + j]).sum::<f64>());
            let gradient: [f64; 2] =
                core::array::from_fn(|j| (0..3).map(|i| g.gradient[i][j] * w[i]).sum());
            0.5 * g.area_m2
                * (stress[0] * gradient[0] * gradient[0]
                    + 2.0 * stress[2] * gradient[0] * gradient[1]
                    + stress[1] * gradient[1] * gradient[1])
        })
        .sum()
}

#[test]
fn g3_signed_anisotropic_geometric_hessian_and_oblique_supports_follow_world_rotation() {
    let original = mesh();
    let angle: f64 = 0.63;
    let (s, c) = angle.sin_cos();
    let rotated_frame = [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]];
    let rotated = ShellMesh::new(
        original
            .nodes
            .iter()
            .map(|p| rotate(rotated_frame, *p))
            .collect(),
        original.tris.clone(),
    )
    .unwrap();
    let stress = [130.0, -40.0, 25.0];
    let free = [-stress[0] / A, -stress[1] / A, -2.0 * stress[2] / A];
    let q = [0.001, -0.002, 0.0005, 0.0015, -0.0007];
    let mut reference: Option<Vec<[f64; 3]>> = None;
    for (mesh, frame) in [(&original, IDENTITY), (&rotated, rotated_frame)] {
        let model = model(mesh, frame, true, false);
        let prepared = model.prepare_thermal(&uniform(free, 4), || false).unwrap();
        let tangent = prepared
            .membrane_prestress_at(&[0.0; 30], limits(), || false)
            .unwrap();
        let force = action(tangent.geometric_stiffness(), &q);
        let energy =
            |v: &[f64]| membrane_energy(mesh, frame[2], &tangent.expand(v).unwrap(), stress);
        near(
            0.5 * q.iter().zip(&force).map(|(a, b)| a * b).sum::<f64>(),
            energy(&q),
            4e-14,
            1e-17,
        );
        for column in 0..q.len() {
            let (mut plus, mut minus) = (q, q);
            let h = 1e-6;
            plus[column] += h;
            minus[column] -= h;
            near(
                (energy(&plus) - energy(&minus)) / (2.0 * h),
                force[column],
                2e-10,
                2e-12,
            );
            for row in 0..q.len() {
                assert_eq!(
                    entry(tangent.geometric_stiffness(), row, column).to_bits(),
                    entry(tangent.geometric_stiffness(), column, row).to_bits()
                );
            }
        }
        let load: Vec<_> = (0..5)
            .flat_map(|node| {
                let f = rotate(frame, [0.0, 0.0, node as f64 + 1.0]);
                [f[0], f[1], f[2], 0.0, 0.0, 0.0]
            })
            .collect();
        for (i, value) in tangent.project_load(&load).unwrap().iter().enumerate() {
            near(*value, i as f64 + 1.0, 3e-15, 1e-14);
        }
        let data: Vec<_> = (0..5)
            .flat_map(|i| (0..5).map(move |j| (i, j)))
            .map(|(i, j)| {
                [
                    entry(tangent.geometric_stiffness(), i, j),
                    entry(&tangent.pencil().k, i, j),
                    entry(&tangent.pencil().m, i, j),
                ]
            })
            .collect();
        if let Some(old) = &reference {
            for (a, b) in data.iter().zip(old) {
                for i in 0..3 {
                    near(a[i], b[i], 2e-12, 1e-10);
                }
            }
        } else {
            reference = Some(data);
        }
    }
}

#[test]
fn g0_prestress_refuses_nonplanar_directors_incompatible_supports_and_large_internal_thermal_rotation()
 {
    let mesh = mesh();
    // Keep the support-frame roundtrip tolerance compatible with the native
    // frame admission, including a state expanded by that exact native map.
    let mut rounded_frame = IDENTITY;
    rounded_frame[2][2] += 4e-13;
    let rounded = model(&mesh, rounded_frame, true, true);
    rounded
        .prepare_thermal(&uniform([0.0; 3], 4), || false)
        .unwrap()
        .membrane_prestress_at(&rounded.expand(&[0.0001]).unwrap(), limits(), || false)
        .unwrap();
    let model = model(&mesh, IDENTITY, true, true);
    let fields = uniform([0.0001; 3], 4);
    let prepared = model.prepare_thermal(&fields, || false).unwrap();
    let mut incompatible = [0.0; 30];
    incompatible[0] = 0.001;
    assert!(
        prepared
            .membrane_prestress_at(&incompatible, limits(), || false)
            .is_err()
    );
    assert!(matches!(
        prepared.membrane_prestress_at(&[0.0; 30], limits(), || true),
        Err(PlateError::Cancelled)
    ));
    let mut too_large = [0.0; 30];
    too_large[6 * 4 + 2] = 0.1;
    assert!(
        prepared
            .membrane_prestress_at(&too_large, limits(), || false)
            .is_err()
    );
    let mut normals = vec![[0.0, 0.0, 1.0]; 5];
    normals[4] = [0.01, 0.0, (1.0_f64 - 0.0001).sqrt()];
    let directors = NodalDirectors::new(&mesh, normals, 0.9).unwrap();
    let variant = assemble_mitc3p(
        &mesh,
        &directors,
        &[section(); 4],
        &[],
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: 5,
            max_triangles: 4,
        },
        1e-3,
        || false,
    )
    .unwrap();
    assert!(
        variant
            .prepare_thermal(&fields, || false)
            .unwrap()
            .membrane_prestress_at(&[0.0; 30], limits(), || false)
            .is_err()
    );
    // A spatial thermal moment excites the original director bubble even when
    // every physical nodal coordinate is zero. Its rotation must be admitted.
    let gradients: Vec<_> = mesh
        .tris
        .iter()
        .map(|tri| {
            FacetThermalStrain::AtNodes(tri.map(|node| {
                ShellThermalStrain::try_new([0.0; 3], [0.005 * mesh.nodes[node][0], 0.0, 0.0])
                    .unwrap()
            }))
        })
        .collect();
    let thermal = model.prepare_thermal(&gradients, || false).unwrap();
    let mut seen = false;
    for element in 0..4 {
        let field = thermal.strain_field(element, &[0.0; 18], || false).unwrap();
        seen |= field
            .reduced_elastic_layer_bound([0.0; 2], [[0.0; 3]; 2], &mut || false)
            .unwrap()
            > 1e-10;
    }
    assert!(seen);
    let strict = MembranePrestressLimits {
        maximum_rotation_rad: 1e-14,
        ..limits()
    };
    assert!(
        thermal
            .membrane_prestress_at(&[0.0; 30], strict, || false)
            .is_err()
    );
    assert_eq!(model.pencil().free, 1);
}
