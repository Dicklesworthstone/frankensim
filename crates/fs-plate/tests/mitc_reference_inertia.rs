//! Native MITC inertia with spatially graded original layer moments.
//! Oracles integrate barycentric monomials independently of element quadrature.

use fs_material::sheet::SheetMatrix;
use fs_plate::PlateError;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::thermal::transient::EnrichedMitcShell;
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellInertia,
    ShearShellSection, ShellConstraint, ShellMass, assemble_mitc3p,
    assemble_mitc3p_with_reference_inertia,
};
use fs_sparse::Csr;
use std::sync::Arc;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn rotate(a: [f64; 3]) -> [f64; 3] {
    [[0.36, -0.48, 0.8], [0.8, 0.6, 0.0], [-0.48, 0.64, 0.6]].map(|r| dot(r, a))
}
fn section(scale: f64) -> ShearShellSection {
    ShearShellSection::new_coupled(
        SheetMatrix::orthotropic(1e6 * scale, 1.4e6 * scale, 0.2, 4e5 * scale).unwrap(),
        [0.0; 9],
        SheetMatrix::orthotropic(400.0 * scale, 500.0 * scale, 0.2, 160.0 * scale).unwrap(),
        [2e5 * scale, 0.0, 2e5 * scale],
        1.0,
        2.0,
        0.01,
        0.002,
        0.0,
    )
    .unwrap()
}
fn triangle() -> (ShellMesh, NodalDirectors) {
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.5).unwrap();
    (mesh, directors)
}
fn budget(mesh: &ShellMesh) -> MitcBudget {
    MitcBudget {
        max_nodes: mesh.nodes.len(),
        max_triangles: mesh.tris.len(),
    }
}
fn model(
    mesh: &ShellMesh,
    directors: &NodalDirectors,
    sections: &[ShearShellSection],
    moments: &[[ShearShellInertia; 3]],
    constraints: &[ShellConstraint],
    mass: ShellMass,
) -> MitcShellModel {
    assemble_mitc3p_with_reference_inertia(
        mesh,
        directors,
        sections,
        moments,
        constraints,
        mass,
        budget(mesh),
        1e-3,
        || false,
    )
    .unwrap()
}
fn entry(m: &Csr, row: usize, col: usize) -> f64 {
    let (cols, values) = m.row(row);
    cols.iter()
        .position(|&c| c == col)
        .map_or(0.0, |i| values[i])
}
fn same_bits(a: &Csr, b: &Csr, n: usize) {
    for row in 0..n {
        let (ac, av) = a.row(row);
        let (bc, bv) = b.row(row);
        assert_eq!(ac, bc);
        assert_eq!(
            av.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            bv.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
    }
}
fn close(observed: f64, expected: f64, scale: f64) {
    assert!(
        (observed - expected).abs() <= 3e-11 * scale,
        "observed={observed:e}, expected={expected:e}, scale={scale:e}"
    );
}
fn rigid(mesh: &ShellMesh, axis: usize) -> Vec<f64> {
    let mut direction = [0.0; 3];
    direction[axis % 3] = 1.0;
    mesh.nodes
        .iter()
        .flat_map(|&p| {
            let (v, w) = if axis < 3 {
                (direction, [0.0; 3])
            } else {
                (cross(direction, p), direction)
            };
            [v[0], v[1], v[2], w[0], w[1], w[2]]
        })
        .collect()
}
fn work(model: &MitcShellModel, a: &[f64], b: &[f64]) -> f64 {
    let (a, b) = (
        model.project_load(a).unwrap(),
        model.project_load(b).unwrap(),
    );
    let mut mb = vec![0.0; b.len()];
    model.pencil().m.spmv(&b, &mut mb);
    a.iter().zip(mb).map(|(a, b)| a * b).sum()
}

#[test]
fn rigid_mass_centroid_and_rotational_inertia_retain_every_corner_moment() {
    for curved in [false, true] {
        for turned in [false, true] {
            for sign in [-1.0, 1.0] {
                let (mut points, tris) = if curved {
                    (
                        vec![
                            [-0.12, -0.07, 0.025],
                            [0.36, -0.04, -0.008],
                            [0.29, 0.31, 0.041],
                            [-0.08, 0.28, 0.018],
                            [0.10, 0.12, 0.09],
                        ],
                        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
                    )
                } else {
                    (
                        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                        vec![[0, 1, 2]],
                    )
                };
                let mut normals: Vec<_> = points
                    .iter()
                    .map(|p| {
                        let d = if curved {
                            [-0.3 * p[0] + 0.1, -0.25 * p[1] - 0.05, 1.0]
                        } else {
                            [0.0, 0.0, 1.0]
                        };
                        d.map(|v| v / dot(d, d).sqrt())
                    })
                    .collect();
                if turned {
                    for p in &mut points {
                        *p = rotate(*p);
                        for (v, shift) in p.iter_mut().zip([0.7, -0.2, 0.4]) {
                            *v += shift;
                        }
                    }
                    for d in &mut normals {
                        *d = rotate(*d);
                    }
                }
                let mesh = ShellMesh::new(points, tris).unwrap();
                let directors = NodalDirectors::new(&mesh, normals.clone(), 0.5).unwrap();
                let sections = vec![section(1.0); mesh.tris.len()];
                let moments: Vec<_> = mesh
                    .tris
                    .iter()
                    .enumerate()
                    .map(|(e, _)| {
                        core::array::from_fn(|i| {
                            let i0 = (1.0 + e as f64) * [1.0, 3.0, 7.0][i];
                            let z = sign * [0.01, -0.02, 0.03][i];
                            ShearShellInertia::new(i0, i0 * z, i0 * (0.0004 + z * z)).unwrap()
                        })
                    })
                    .collect();
                let assembled = model(
                    &mesh,
                    &directors,
                    &sections,
                    &moments,
                    &[],
                    ShellMass::Consistent,
                );
                let mut exact = [[0.0; 6]; 6];
                for (e, tri) in mesh.tris.iter().enumerate() {
                    let edge = |i: usize| {
                        core::array::from_fn(|k| mesh.nodes[tri[i]][k] - mesh.nodes[tri[0]][k])
                    };
                    let normal = cross(edge(1), edge(2));
                    let area = 0.5 * dot(normal, normal).sqrt();
                    let fields: [([[f64; 3]; 3], [[f64; 3]; 3]); 6] = core::array::from_fn(|a| {
                        let mut axis = [0.0; 3];
                        axis[a % 3] = 1.0;
                        if a < 3 {
                            ([axis; 3], [[0.0; 3]; 3])
                        } else {
                            (
                                tri.map(|i| cross(axis, mesh.nodes[i])),
                                tri.map(|i| cross(axis, normals[i])),
                            )
                        }
                    });
                    for a in 0..6 {
                        for b in 0..6 {
                            for i in 0..3 {
                                for j in 0..3 {
                                    for k in 0..3 {
                                        // integral lambda_i lambda_j lambda_k = 2A*prod(alpha!)/5!.
                                        let factor = if i == j && j == k {
                                            6.0
                                        } else if i == j || i == k || j == k {
                                            2.0
                                        } else {
                                            1.0
                                        };
                                        let moment = moments[e][k];
                                        exact[a][b] += area / 60.0
                                            * factor
                                            * (moment.areal_mass_kg_m2()
                                                * dot(fields[a].0[i], fields[b].0[j])
                                                + moment.mass_first_moment_kg_m()
                                                    * (dot(fields[a].0[i], fields[b].1[j])
                                                        + dot(fields[a].1[i], fields[b].0[j]))
                                                + moment.rotary_inertia_kg()
                                                    * dot(fields[a].1[i], fields[b].1[j]));
                                    }
                                }
                            }
                        }
                    }
                }
                for a in 0..6 {
                    for b in 0..6 {
                        close(
                            work(&assembled, &rigid(&mesh, a), &rigid(&mesh, b)),
                            exact[a][b],
                            (exact[a][a] * exact[b][b]).sqrt(),
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn constant_corner_moments_preserve_both_native_mass_paths_bitwise() {
    let (mesh, directors) = triangle();
    let sections = [section(1.0)];
    let moments = [[ShearShellInertia::from_section(&sections[0]); 3]];
    for mass in [ShellMass::Consistent, ShellMass::Lumped] {
        let original = assemble_mitc3p(
            &mesh,
            &directors,
            &sections,
            &[],
            mass,
            budget(&mesh),
            1e-3,
            || false,
        )
        .unwrap();
        let explicit = model(&mesh, &directors, &sections, &moments, &[], mass);
        same_bits(&original.pencil().k, &explicit.pencil().k, 18);
        same_bits(&original.pencil().m, &explicit.pencil().m, 18);
        assert_eq!(
            original.reference_inertia_moments(0),
            explicit.reference_inertia_moments(0)
        );
        assert!(explicit.reference_inertia_moments(1).is_none());
    }
}

fn scaled_moments(mass: [f64; 3]) -> [[ShearShellInertia; 3]; 1] {
    [mass.map(|m| ShearShellInertia::new(m, 0.01 * m, 0.001 * m).unwrap())]
}
fn single_dof() -> Vec<ShellConstraint> {
    (0..3)
        .map(|node| {
            let mut mask = [true; 6];
            if node == 0 {
                mask[0] = false;
            }
            ShellConstraint {
                node,
                mask,
                frame: ConstraintFrame::Global,
            }
        })
        .collect()
}

#[test]
fn same_total_mass_with_different_placement_changes_native_modal_inertia() {
    let (mesh, directors) = triangle();
    let sections = [section(1.0)];
    let profiles = [
        scaled_moments([1.0, 2.0, 5.0]),
        scaled_moments([5.0, 2.0, 1.0]),
    ];
    let mut eigenvalues = Vec::new();
    for (field, exact_mass) in profiles.iter().zip([1.0 / 6.0, 0.3]) {
        let m = model(
            &mesh,
            &directors,
            &sections,
            field,
            &single_dof(),
            ShellMass::Consistent,
        );
        assert_eq!(m.pencil().free, 1);
        // For the free in-plane N0 coordinate, bubble stiffness coupling is zero.
        // integral I0*N0² = A*(3*I0_0+I0_1+I0_2)/30, while total mass is 4/3.
        close(entry(&m.pencil().m, 0, 0), exact_mass, exact_mass);
        eigenvalues.push(entry(&m.pencil().k, 0, 0) / entry(&m.pencil().m, 0, 0));
    }
    close(eigenvalues[0] / eigenvalues[1], 1.8, 1.8);
}

#[test]
fn retained_director_dynamics_use_the_original_profile_after_material_refresh() {
    let (mesh, directors) = triangle();
    let field = scaled_moments([1.0, 2.0, 5.0]);
    let original = model(
        &mesh,
        &directors,
        &[section(1.0)],
        &field,
        &single_dof(),
        ShellMass::Consistent,
    );
    let enriched = EnrichedMitcShell::try_new(original, || false).unwrap();
    assert_eq!(enriched.dimension(), 3);
    close(entry(enriched.mass(), 0, 0), 1.0 / 6.0, 1.0 / 6.0);
    // Independent degree-seven integral: b=27*N0*N1*N2,
    // integral I2*b² = (27*A/280)*sum(corner I2).
    let exact_bubble = 27.0 * 0.5 / 280.0 * 0.008;
    for i in 1..3 {
        close(entry(enriched.mass(), i, i), exact_bubble, exact_bubble);
    }
    // integral I1*N0*b = A*(9*I1_0/140 + 3*(I1_1+I1_2)/70).
    // Squared norm is independent of the element's two orthonormal bubble axes.
    let cross_moment = 0.5 * (9.0 * 0.01 / 140.0 + 3.0 * 0.07 / 70.0);
    let observed =
        (entry(enriched.mass(), 0, 1).powi(2) + entry(enriched.mass(), 0, 2).powi(2)).sqrt();
    close(observed, cross_moment, cross_moment);
    let current = model(
        &mesh,
        &directors,
        &[section(2.0)],
        &field,
        &single_dof(),
        ShellMass::Consistent,
    );
    let refreshed = enriched
        .with_material_model(Arc::new(current), || false)
        .unwrap();
    same_bits(enriched.mass(), refreshed.mass(), enriched.dimension());
    let relocated = model(
        &mesh,
        &directors,
        &[section(2.0)],
        &scaled_moments([5.0, 2.0, 1.0]),
        &single_dof(),
        ShellMass::Consistent,
    );
    assert!(
        enriched
            .with_material_model(Arc::new(relocated), || false)
            .is_err()
    );
}

#[test]
fn lumped_mass_uses_exact_positive_corner_partitions() {
    let (mesh, directors) = triangle();
    let field = scaled_moments([1.0, 2.0, 5.0]);
    let m = model(
        &mesh,
        &directors,
        &[section(1.0)],
        &field,
        &[],
        ShellMass::Lumped,
    );
    for (node, value) in [1.0, 2.0, 5.0].into_iter().enumerate() {
        let i0 = 0.5 / 12.0 * (8.0 + value);
        for axis in 0..3 {
            close(
                entry(&m.pencil().m, 6 * node + axis, 6 * node + axis),
                i0,
                i0,
            );
            close(
                entry(&m.pencil().m, 6 * node + 3 + axis, 6 * node + 3 + axis),
                0.001 * i0,
                0.001 * i0,
            );
        }
        close(
            entry(&m.pencil().m, 6 * node, 6 * node + 4),
            0.01 * i0,
            0.01 * i0,
        );
        close(
            entry(&m.pencil().m, 6 * node + 1, 6 * node + 3),
            -0.01 * i0,
            0.01 * i0,
        );
    }
}

#[test]
fn unrealizable_inventory_wrong_extent_and_cancellation_refuse() {
    for [i0, i1, i2] in [
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
        [2.0, 2.0, 2.0],
        [2.0, -2.0, 2.0],
        [2.0, f64::from_bits(2.0_f64.to_bits() + 1), 2.0],
        [2.0, -f64::from_bits(2.0_f64.to_bits() + 1), 2.0],
        [f64::MAX, f64::MAX, f64::MAX],
        [f64::from_bits(1); 3],
        [f64::MIN_POSITIVE, 3.0, f64::MAX],
        [1.0, -2.0, 1.0],
        [f64::NAN, 0.0, 1.0],
        [1.0, 0.0, f64::INFINITY],
    ] {
        assert!(ShearShellInertia::new(i0, i1, i2).is_err());
    }
    // Exact positive blocks stay admitted even when floating products would
    // overflow/underflow, or the determinant lies below a heuristic margin.
    for [i0, i1, i2] in [
        [f64::from_bits(1), 0.0, f64::from_bits(1)],
        [f64::MAX, f64::MAX / 2.0, f64::MAX],
        [f64::MIN_POSITIVE, 1.0, f64::MAX],
        [1.0, 1.0 - f64::EPSILON, 1.0],
        [2.0, f64::from_bits(2.0_f64.to_bits() - 1), 2.0],
        [2.0, -f64::from_bits(2.0_f64.to_bits() - 1), 2.0],
    ] {
        assert!(ShearShellInertia::new(i0, i1, i2).is_ok());
    }
    let (mesh, directors) = triangle();
    assert!(matches!(
        assemble_mitc3p_with_reference_inertia(
            &mesh,
            &directors,
            &[section(1.0)],
            &[],
            &[],
            ShellMass::Consistent,
            budget(&mesh),
            1e-3,
            || false
        ),
        Err(PlateError::SectionCount {
            expected: 1,
            actual: 0
        })
    ));
    assert!(matches!(
        assemble_mitc3p_with_reference_inertia(
            &mesh,
            &directors,
            &[section(1.0)],
            &scaled_moments([1.0, 2.0, 5.0]),
            &[],
            ShellMass::Consistent,
            budget(&mesh),
            1e-3,
            || true
        ),
        Err(PlateError::Cancelled)
    ));
}
