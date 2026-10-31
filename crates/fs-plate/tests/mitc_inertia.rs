//! Rigid spatial inertia of the original condensed MITC reference-area/director
//! field. Section moments remain explicit; no curved-volume Jacobian is inferred.

use fs_material::sheet::SheetMatrix;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellSection,
    ShellConstraint, ShellMass, assemble_mitc3p,
};
use fs_sparse::{DirectOrdering, LdltOptions, SymbolicLdlt};

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
fn facet(mesh: &ShellMesh, tri: [usize; 3]) -> (f64, [f64; 3]) {
    let edge = |i: usize| core::array::from_fn(|k| mesh.nodes[i][k] - mesh.nodes[tri[0]][k]);
    let n = cross(edge(tri[1]), edge(tri[2]));
    let length = dot(n, n).sqrt();
    (0.5 * length, n.map(|x| x / length))
}
fn section(index: usize, sign: f64) -> ShearShellSection {
    let a = SheetMatrix::orthotropic(1e6, 1.4e6, 0.2, 4e5).unwrap();
    let d = SheetMatrix::orthotropic(400., 500., 0.2, 160.).unwrap();
    let i0 = 2.0 + 0.1 * index as f64;
    ShearShellSection::new_coupled(
        a,
        [0.; 9],
        d,
        [2e5, 0., 2e5],
        1.,
        i0,
        sign * 0.0045 * i0,
        0.00008 * i0,
        0.,
    )
    .unwrap()
}
fn fixture(curved: bool, turned: bool) -> (ShellMesh, Vec<[f64; 3]>) {
    let (mut points, triangles) = if curved {
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
            vec![[0., 0., 0.], [0.35, 0., 0.], [0.04, 0.27, 0.]],
            vec![[0, 1, 2]],
        )
    };
    let mut directors: Vec<_> = points
        .iter()
        .map(|p| {
            let n = if curved {
                [-0.3 * p[0] + 0.1, -0.25 * p[1] - 0.05, 1.]
            } else {
                [0., 0., 1.]
            };
            let norm = dot(n, n).sqrt();
            n.map(|x| x / norm)
        })
        .collect();
    if turned {
        for p in &mut points {
            *p = rotate(*p);
            for (x, shift) in p.iter_mut().zip([0.7, -0.2, 0.4]) {
                *x += shift;
            }
        }
        for n in &mut directors {
            *n = rotate(*n);
        }
    }
    (ShellMesh::new(points, triangles).unwrap(), directors)
}
fn assemble(
    mesh: &ShellMesh,
    directors: &[[f64; 3]],
    sections: &[ShearShellSection],
    constraints: &[ShellConstraint],
    mass: ShellMass,
) -> MitcShellModel {
    let directors = NodalDirectors::new(mesh, directors.to_vec(), 0.5).unwrap();
    assemble_mitc3p(
        mesh,
        &directors,
        sections,
        constraints,
        mass,
        MitcBudget {
            max_nodes: mesh.nodes.len(),
            max_triangles: mesh.tris.len(),
        },
        1e-3,
        || false,
    )
    .unwrap()
}
fn rigid(mesh: &ShellMesh, axis: usize) -> Vec<f64> {
    let mut direction = [0.; 3];
    direction[axis % 3] = 1.;
    mesh.nodes
        .iter()
        .flat_map(|&p| {
            let (v, w) = if axis < 3 {
                (direction, [0.; 3])
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
    let mut mb = vec![0.; b.len()];
    model.pencil().m.spmv(&b, &mut mb);
    a.iter().zip(mb).map(|(a, b)| a * b).sum()
}

// Independent barycentric volume-moment integration: for P1 vector fields,
// integral a dot b = area/12 * [(sum a_i) dot (sum b_i) + sum a_i dot b_i].
// Rigid director motion is omega cross d; all three physical I0/I1/I2 terms
// are retained, including translation/rotation coupling about a remote origin.
fn rigid_moments(
    mesh: &ShellMesh,
    directors: &[[f64; 3]],
    sections: &[ShearShellSection],
) -> [[f64; 6]; 6] {
    let mut out = [[0.; 6]; 6];
    for (e, &tri) in mesh.tris.iter().enumerate() {
        let area = facet(mesh, tri).0;
        let s = sections[e];
        let fields: [([[f64; 3]; 3], [[f64; 3]; 3]); 6] = core::array::from_fn(|a| {
            let mut axis = [0.; 3];
            axis[a % 3] = 1.;
            if a < 3 {
                ([axis; 3], [[0.; 3]; 3])
            } else {
                (
                    tri.map(|i| cross(axis, mesh.nodes[i])),
                    tri.map(|i| cross(axis, directors[i])),
                )
            }
        });
        for a in 0..6 {
            for b in 0..6 {
                for i in 0..3 {
                    for j in 0..3 {
                        let gram = area / 12. * if i == j { 2. } else { 1. };
                        out[a][b] += gram
                            * (s.areal_mass_kg_m2() * dot(fields[a].0[i], fields[b].0[j])
                                + s.mass_first_moment_kg_m()
                                    * (dot(fields[a].0[i], fields[b].1[j])
                                        + dot(fields[a].1[i], fields[b].0[j]))
                                + s.rotary_inertia_kg() * dot(fields[a].1[i], fields[b].1[j]));
                    }
                }
            }
        }
    }
    out
}

#[test]
fn full_rigid_spatial_inertia_matches_section_moments_on_flat_and_curved_fans() {
    for curved in [false, true] {
        for turned in [false, true] {
            for sign in [-1., 1.] {
                let (mesh, directors) = fixture(curved, turned);
                let sections: Vec<_> = (0..mesh.tris.len()).map(|e| section(e, sign)).collect();
                let model = assemble(&mesh, &directors, &sections, &[], ShellMass::Consistent);
                let expected = rigid_moments(&mesh, &directors, &sections);
                let fields: Vec<_> = (0..6).map(|axis| rigid(&mesh, axis)).collect();
                for a in 0..6 {
                    for b in 0..6 {
                        let observed = work(&model, &fields[a], &fields[b]);
                        let scale = (expected[a][a] * expected[b][b]).sqrt();
                        assert!(
                            (observed - expected[a][b]).abs() <= 3e-11 * scale,
                            "curved={curved} turned={turned} sign={sign} [{a},{b}]: {observed:e} != {:e}",
                            expected[a][b]
                        );
                    }
                }
                // The previous absolute drilling term added area*I2 to normal spin,
                // even though a flat plate's material director does not move at all.
                if !curved && !turned {
                    let old_excess = facet(&mesh, mesh.tris[0]).0 * sections[0].rotary_inertia_kg();
                    assert!(old_excess > 1e-4 * expected[5][5]);
                }
            }
        }
    }
}

#[test]
fn relative_gauge_is_positive_and_orthogonal_to_every_rigid_velocity() {
    for turned in [false, true] {
        let (mesh, directors) = fixture(true, turned);
        let sections: Vec<_> = (0..mesh.tris.len()).map(|e| section(e, 1.)).collect();
        let model = assemble(&mesh, &directors, &sections, &[], ShellMass::Consistent);
        let factors: Vec<_> = (0..mesh.nodes.len())
            .map(|i| 0.3 + 0.2 * i as f64)
            .collect();
        let gauge: Vec<_> = directors
            .iter()
            .zip(&factors)
            .flat_map(|(d, b)| [0., 0., 0., b * d[0], b * d[1], b * d[2]])
            .collect();
        let mut areas = vec![0.; mesh.nodes.len()];
        let mut normals = vec![[0.; 3]; mesh.nodes.len()];
        let mut inertia = vec![0.; mesh.nodes.len()];
        for (e, &tri) in mesh.tris.iter().enumerate() {
            let (area, normal) = facet(&mesh, tri);
            for i in tri {
                areas[i] += area / 3.;
                inertia[i] += area / 3. * sections[e].rotary_inertia_kg();
                for c in 0..3 {
                    normals[i][c] += area / 3. * normal[c];
                }
            }
        }
        let expected: f64 = (0..mesh.nodes.len())
            .map(|i| inertia[i] * (factors[i] * dot(normals[i], directors[i]) / areas[i]).powi(2))
            .sum();
        let actual = work(&model, &gauge, &gauge);
        assert!((actual / expected - 1.).abs() < 3e-11);
        for axis in 0..6 {
            let r = rigid(&mesh, axis);
            let scale = (work(&model, &r, &r) * actual).sqrt();
            assert!(work(&model, &r, &gauge).abs() < 3e-11 * scale);
        }
        for constraints in [
            vec![],
            vec![ShellConstraint {
                node: 0,
                mask: [true; 6],
                frame: ConstraintFrame::Global,
            }],
        ] {
            let m = assemble(
                &mesh,
                &directors,
                &sections,
                &constraints,
                ShellMass::Consistent,
            );
            let factor = SymbolicLdlt::analyze(&m.pencil().m, DirectOrdering::Amd)
                .unwrap()
                .factor(&m.pencil().m, &LdltOptions::default())
                .unwrap();
            assert_eq!(factor.inertia().positive, m.pencil().free);
            assert_eq!(factor.inertia().negative, 0);
        }
    }
}

#[test]
fn explicit_lumped_inertia_keeps_original_nodal_matrix_and_stiffness() {
    let (mesh, directors) = fixture(false, false);
    let sections = [section(0, 1.)];
    let lumped = assemble(&mesh, &directors, &sections, &[], ShellMass::Lumped);
    let consistent = assemble(&mesh, &directors, &sections, &[], ShellMass::Consistent);
    assert_eq!(lumped.pencil().k, consistent.pencil().k);
    let mut expected = [[0.; 18]; 18];
    let weight = facet(&mesh, mesh.tris[0]).0 / 3.;
    for i in 0..3 {
        for c in 0..6 {
            expected[6 * i + c][6 * i + c] = weight
                * if c < 3 {
                    sections[0].areal_mass_kg_m2()
                } else {
                    sections[0].rotary_inertia_kg()
                };
        }
        for c in 0..3 {
            let mut axis = [0.; 3];
            axis[c] = 1.;
            let motion = cross(axis, directors[i]);
            for d in 0..3 {
                let value = weight * sections[0].mass_first_moment_kg_m() * motion[d];
                expected[6 * i + d][6 * i + 3 + c] = value;
                expected[6 * i + 3 + c][6 * i + d] = value;
            }
        }
    }
    for i in 0..18 {
        for j in 0..18 {
            assert_eq!(lumped.pencil().m.get(i, j), expected[i][j]);
        }
    }
}
