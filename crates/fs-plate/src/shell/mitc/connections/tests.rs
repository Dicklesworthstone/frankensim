use super::*;
use crate::shell::{ShellMesh, mitc::{assemble_mitc3p, ConstraintFrame, MitcBudget,
    NodalDirectors, ShearShellSection, ShellConstraint, ShellMass, solve_static}};

fn model() -> MitcShellModel {
    let mesh = ShellMesh::new(vec![[0.0, 0.0, 0.0], [0.1, 0.0, 0.0],
        [0.0, 0.1, 0.0]], vec![[0, 1, 2]]).unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.99).unwrap();
    let section = ShearShellSection::from_orthotropic_solid(1e9, 4e8, 0.2, 2e8,
        1e8, 0.6e8, 0.003, 700.0, 0.17, 5.0 / 6.0).unwrap();
    let s = fs_math::det::sqrt(0.5);
    let supports = [ShellConstraint { node: 0, mask: [true; 6], frame: ConstraintFrame::Global },
        ShellConstraint { node: 1, mask: [true; 6], frame: ConstraintFrame::Global },
        ShellConstraint { node: 2, mask: [false; 6],
            frame: ConstraintFrame::Local([[1.0, 0.0, 0.0], [0.0, s, s], [0.0, -s, s]]) }];
    assemble_mitc3p(&mesh, &directors, &[section], &supports, ShellMass::Consistent,
        MitcBudget { max_nodes: 3, max_triangles: 1 }, 1e-3, || false).unwrap()
}
fn budget() -> ShellConnectionBudget {
    ShellConnectionBudget { max_connections: 2, max_terms: 4, max_work: 100_000 }
}
fn link() -> ShellConnection {
    ShellConnection { terms: vec![ShellConnectionTerm { node: 2,
        translation: [0.0, 0.0, 1.0], rotation_m: [0.0, -0.03, 0.0] }],
        stiffness_n_m: 3000.0, damping_n_s_m: 0.7 }
}

#[test]
fn oblique_support_and_offset_attachment_match_independent_static_compliance() {
    let original = model();
    let connection = link();
    let changed = original.with_connections(&[connection.clone()], budget(), || false).unwrap();
    let force = physical_row(&connection, 3).unwrap();
    let displacement = solve_static(&original, &force).unwrap();
    let after = solve_static(&changed, &force).unwrap();
    let compliance: f64 = force.iter().zip(&displacement).map(|(a, b)| a * b).sum();
    let actual: f64 = force.iter().zip(&after).map(|(a, b)| a * b).sum();
    // Independently eliminate the scalar reaction F=-k*u from u=C*(1+F).
    let expected = compliance / (1.0 + connection.stiffness_n_m * compliance);
    assert!(compliance > 0.0 && actual > 0.0 && actual < compliance);
    assert!((actual / expected - 1.0).abs() < 1e-8);
    for i in 0..original.pencil.free {
        assert_eq!(original.pencil.m.row(i), changed.pencil.m.row(i));
    }
    let free = vec![0.03, -0.02, 0.01, -0.4, 0.2, 0.1];
    let physical = changed.expand(&free).unwrap();
    let projected = changed.project_connections(&[connection.clone()], &[physical.clone()],
        budget(), || false).unwrap()[0][0];
    let direct = physical[14] - 0.03 * physical[16];
    assert!((projected - direct).abs() < 1e-14);
    let mut before = vec![0.0; free.len()];
    let mut after = before.clone();
    original.pencil.k.spmv(&free, &mut before);
    changed.pencil.k.spmv(&free, &mut after);
    let added: f64 = free.iter().zip(before.iter().zip(&after)).map(|(q, (a, b))| q * (b - a)).sum();
    let expected = connection.stiffness_n_m * projected * projected;
    assert!((added - expected).abs() < 1e-9 * expected.max(1.0));
}

#[test]
fn coincident_offset_sides_preserve_common_rigid_motion_and_viscous_work() {
    let model = model();
    // Nodes 1 and 2 meet at the same physical point through explicit offsets.
    // Rotation is essential: a translation-only tie would penalize rigid motion.
    let link = ShellConnection { terms: vec![
        ShellConnectionTerm { node: 1, translation: [0.0, 0.0, 1.0],
            rotation_m: [0.05, 0.05, 0.0] },
        ShellConnectionTerm { node: 2, translation: [0.0, 0.0, -1.0],
            rotation_m: [0.05, 0.05, 0.0] },
    ], stiffness_n_m: 1000.0, damping_n_s_m: 2.0 };
    let positions = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.1, 0.0]];
    let translation = [0.2, 0.3, -0.1];
    let theta = [0.4, -0.7, 0.1];
    let mut rigid = vec![0.0; 18];
    for (i, p) in positions.iter().enumerate() {
        let rotation = [theta[1] * p[2] - theta[2] * p[1],
            theta[2] * p[0] - theta[0] * p[2], theta[0] * p[1] - theta[1] * p[0]];
        for c in 0..3 { rigid[6 * i + c] = translation[c] + rotation[c]; rigid[6 * i + 3 + c] = theta[c]; }
    }
    let row = model.project_connections(&[link.clone()], &[rigid], budget(), || false).unwrap();
    assert!(row[0][0].abs() < 1e-14);
    let mut moving = vec![0.0; 18];
    moving[8] = 0.2; moving[14] = -0.1; moving[9] = 0.4;
    let velocity = model.project_connections(&[link.clone()], &[moving.clone()],
        budget(), || false).unwrap()[0][0];
    let force = physical_row(&link, 3).unwrap();
    let work: f64 = force.iter().zip(moving).map(|(b, v)| -link.damping_n_s_m * velocity * b * v).sum();
    assert!((work + link.damping_n_s_m * velocity * velocity).abs() < 1e-14);
    assert!(work < 0.0);
}

#[test]
fn connection_refusals_do_not_publish_partial_stiffness() {
    let model = model();
    let mut connection = link();
    let original = model.pencil.k.row(0).1.to_vec();
    let mut too_small = budget(); too_small.max_work = 1;
    assert!(model.with_connections(&[connection.clone()], too_small, || false).is_err());
    assert!(model.with_connections(&[connection.clone()], budget(), || true).is_err());
    connection.terms[0].node = 0;
    assert!(model.with_connections(&[connection.clone()], budget(), || false).is_err());
    connection.terms[0].node = 2; connection.damping_n_s_m = -1.0;
    assert!(model.with_connections(&[connection], budget(), || false).is_err());
    assert_eq!(original, model.pencil.k.row(0).1);
}
