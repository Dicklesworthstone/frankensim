//! Physical CST/DKT material loss, independently checked in strain coordinates.

use fs_plate::shell::viscosity::{
    ShellViscosity, ShellViscosityBudget, assemble_constitutive_operator, project_viscosity,
    projection_terms,
};
use fs_plate::shell::{ShellMesh, ShellModel, ShellSupport, assemble_shell};
use fs_plate::{PlateSection, dkt_stiffness};

fn mesh() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}
fn model(mesh: &ShellMesh) -> ShellModel {
    let section = PlateSection::isotropic(1e7, 0.25, 0.01, 800.0).unwrap();
    assemble_shell(mesh, &section, &[], ShellSupport::Free).unwrap()
}
fn budget() -> ShellViscosityBudget {
    ShellViscosityBudget {
        max_elements: 4,
        max_full_dofs: 36,
        max_dofs: 36,
        max_modes: 18,
        max_work: 20_000_000,
    }
}
fn bilinear(a: &[f64; 9], x: [f64; 3], y: [f64; 3]) -> f64 {
    (0..3)
        .map(|i| (0..3).map(|j| x[i] * a[3 * i + j] * y[j]).sum::<f64>())
        .sum()
}
fn matrix_work(matrix: &fs_sparse::Csr, left: &[f64], right: &[f64]) -> f64 {
    let mut force = vec![0.0; right.len()];
    matrix.spmv(right, &mut force);
    left.iter().zip(force).map(|(x, y)| x * y).sum()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-12 * expected.abs().max(1e-12),
        "{actual:e} != {expected:e}"
    );
}
fn patch(mesh: &ShellMesh, membrane: [f64; 3], kappa: [f64; 3]) -> Vec<f64> {
    mesh.nodes
        .iter()
        .flat_map(|p| {
            let [x, y, _] = *p;
            [
                membrane[0] * x + 0.5 * membrane[2] * y,
                membrane[1] * y + 0.5 * membrane[2] * x,
                -0.5 * kappa[0] * x * x - 0.5 * kappa[1] * y * y - 0.5 * kappa[2] * x * y,
                -kappa[1] * y - 0.5 * kappa[2] * x,
                kappa[0] * x + 0.5 * kappa[2] * y,
                0.0,
            ]
        })
        .collect()
}

#[test]
fn rigid_motion_and_drilling_dissipate_no_material_energy() {
    let mesh = mesh();
    let model = model(&mesh);
    let law = ShellViscosity::new(
        [6.0, 1.0, 0.0, 1.0, 5.0, 0.0, 0.0, 0.0, 2.0],
        [0.8, 0.1, 0.0, 0.1, 0.4, 0.0, 0.0, 0.0, 0.3],
    )
    .unwrap();
    let c =
        assemble_constitutive_operator(&mesh, &[law.clone()], &model, budget(), || false).unwrap();
    let mut vectors = Vec::new();
    for axis in 0..6 {
        let mut v = vec![0.0; 18];
        for (node, p) in mesh.nodes.iter().enumerate() {
            if axis < 3 {
                v[6 * node + axis] = 1.0;
            } else {
                let mut omega = [0.0; 3];
                omega[axis - 3] = 1.0;
                v[6 * node] = omega[1] * p[2] - omega[2] * p[1];
                v[6 * node + 1] = omega[2] * p[0] - omega[0] * p[2];
                v[6 * node + 2] = omega[0] * p[1] - omega[1] * p[0];
                v[6 * node + 3..6 * node + 6].copy_from_slice(&omega);
            }
        }
        vectors.push(v);
    }
    let mut drilling = vec![0.0; 18];
    drilling[5] = 1.0;
    drilling[11] = -3.0;
    drilling[17] = 2.0;
    assert!(
        matrix_work(&model.k, &drilling, &drilling) > 0.0,
        "the elastic numerical drilling gauge must actually be active"
    );
    vectors.push(drilling);
    for vector in &vectors {
        let mut force = vec![0.0; 18];
        c.spmv(vector, &mut force);
        assert!(force.iter().all(|v| v.abs() < 4e-14));
    }
    let projected = project_viscosity(&mesh, &[law], &model, &vectors, budget(), || false).unwrap();
    assert!(projected.iter().all(|v| v.abs() < 1e-28));
    assert_eq!(
        projected[7 * 6 + 6],
        0.0,
        "relative spin must receive exactly zero physical loss"
    );
}

#[test]
fn affine_membrane_and_bending_match_physical_layer_work_with_signed_coupling() {
    let mesh = mesh();
    let model = model(&mesh);
    let a = [8.0, 1.0, 0.0, 1.0, 6.0, 0.0, 0.0, 0.0, 4.0];
    let b = [0.5, 0.1, 0.0, 0.1, 0.2, 0.0, 0.0, 0.0, 0.3];
    let d = [5.0, 0.2, 0.0, 0.2, 3.0, 0.0, 0.0, 0.0, 2.0];
    let strains = [
        ([0.3, -0.2, 0.1], [0.4, 0.2, -0.3]),
        ([-0.1, 0.5, -0.2], [-0.3, 0.1, 0.4]),
    ];
    let vectors: Vec<_> = strains.iter().map(|(e, k)| patch(&mesh, *e, *k)).collect();
    for sign in [1.0, -1.0] {
        let coupling = b.map(|x| sign * x);
        let law = ShellViscosity::new_coupled(a, coupling, d).unwrap();
        let c = assemble_constitutive_operator(&mesh, &[law.clone()], &model, budget(), || false)
            .unwrap();
        let p = project_viscosity(&mesh, &[law], &model, &vectors, budget(), || false).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                let (ei, ki) = strains[i];
                let (ej, kj) = strains[j];
                let expected = 0.5
                    * (bilinear(&a, ei, ej)
                        + bilinear(&coupling, ei, kj)
                        + bilinear(&coupling, ki, ej)
                        + bilinear(&d, ki, kj));
                near(p[2 * i + j], expected);
                near(matrix_work(&c, &vectors[i], &vectors[j]), expected);
                assert_eq!(p[2 * i + j], p[2 * j + i]);
            }
        }
    }
}

#[test]
fn pure_bending_matches_the_original_flat_plate_dkt_operator() {
    let mesh = mesh();
    let model = model(&mesh);
    let d = [4.0, 0.5, 0.2, 0.5, 3.0, 0.1, 0.2, 0.1, 2.0];
    let law = ShellViscosity::new([0.0; 9], d).unwrap();
    let c = assemble_constitutive_operator(&mesh, &[law], &model, budget(), || false).unwrap();
    let (plate, _) = dkt_stiffness(&[0.0, 1.0, 0.0], &[0.0, 0.0, 1.0], &d, 0).unwrap();
    let slots = [2, 4, 3, 8, 10, 9, 14, 16, 15];
    let signs = [1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0];
    let scale = plate.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    for i in 0..9 {
        for j in 0..9 {
            let expected = signs[i] * plate[9 * i + j] * signs[j];
            assert!((c.get(slots[i], slots[j]) - expected).abs() < 2e-13 * scale);
        }
    }
    for i in 0..18 {
        for j in 0..18 {
            if !slots.contains(&i) || !slots.contains(&j) {
                assert_eq!(c.get(i, j), 0.0);
            }
        }
    }
}

#[test]
fn facet_rotation_and_original_support_elimination_preserve_material_work() {
    let mesh = mesh();
    let rotate = |p: [f64; 3]| [0.8 * p[0] + 0.6 * p[2], p[1], -0.6 * p[0] + 0.8 * p[2]];
    let tilted = ShellMesh::new(
        mesh.nodes.iter().copied().map(rotate).collect(),
        mesh.tris.clone(),
    )
    .unwrap();
    let law = ShellViscosity::new(
        [3.0, 0.4, 0.1, 0.4, 2.0, 0.2, 0.1, 0.2, 1.0],
        [0.5, 0.1, 0.0, 0.1, 0.3, 0.0, 0.0, 0.0, 0.2],
    )
    .unwrap();
    let mut vector = patch(&mesh, [0.3, 0.1, -0.2], [0.2, 0.4, -0.1]);
    // A fully clamped node at the origin is exactly compatible with this field.
    vector[..6].fill(0.0);
    let rotated: Vec<_> = vector
        .chunks_exact(3)
        .flat_map(|p| rotate([p[0], p[1], p[2]]))
        .collect();
    let section = PlateSection::isotropic(1e7, 0.25, 0.01, 800.0).unwrap();
    let flat_model = assemble_shell(&mesh, &section, &[0], ShellSupport::Clamped).unwrap();
    let tilt_model = assemble_shell(&tilted, &section, &[0], ShellSupport::Clamped).unwrap();
    let flat = project_viscosity(
        &mesh,
        &[law.clone()],
        &flat_model,
        &[vector[6..].to_vec()],
        budget(),
        || false,
    )
    .unwrap();
    let tilt = project_viscosity(
        &tilted,
        &[law],
        &tilt_model,
        &[rotated[6..].to_vec()],
        budget(),
        || false,
    )
    .unwrap();
    near(tilt[0], flat[0]);
}

#[test]
fn weak_directional_loss_and_exact_nulls_survive_direct_projection() {
    let mesh = mesh();
    let model = model(&mesh);
    let diagonal = 1.0 + 1e-14;
    let law =
        ShellViscosity::new([1.0, 1.0, 0.0, 1.0, diagonal, 0.0, 0.0, 0.0, 0.0], [0.0; 9]).unwrap();
    let amplitudes = [1.0, -2.0, 0.0];
    let vectors: Vec<_> = amplitudes
        .iter()
        .map(|a| patch(&mesh, [-*a, *a, 0.0], [0.0; 3]))
        .collect();
    let p = project_viscosity(&mesh, &[law], &model, &vectors, budget(), || false).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            let expected = 0.5 * (diagonal - 1.0) * amplitudes[i] * amplitudes[j];
            if expected == 0.0 {
                assert_eq!(p[3 * i + j], 0.0);
            } else {
                assert!((p[3 * i + j] / expected - 1.0).abs() < 1e-12);
            }
        }
    }
}

#[test]
fn incomplete_or_active_laws_maps_budgets_and_cancellation_refuse() {
    let mesh = mesh();
    let model = model(&mesh);
    let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let law = ShellViscosity::new(identity, identity).unwrap();
    assert!(ShellViscosity::new(identity.map(|x| -x), identity).is_err());
    assert!(ShellViscosity::new_coupled(identity, identity.map(|x| 2.0 * x), identity).is_err());
    assert!(assemble_constitutive_operator(&mesh, &[], &model, budget(), || false).is_err());
    let mut stale = model.clone();
    stale.dof_map[1] = Some(0);
    assert!(
        assemble_constitutive_operator(&mesh, &[law.clone()], &stale, budget(), || false).is_err()
    );
    let vectors = vec![vec![1.0; 18]];
    let exact = projection_terms(18, 1, 1).unwrap();
    let at_limit = ShellViscosityBudget {
        max_work: exact,
        ..budget()
    };
    assert!(project_viscosity(&mesh, &[law.clone()], &model, &vectors, at_limit, || false).is_ok());
    let too_small = ShellViscosityBudget {
        max_work: exact - 1,
        ..budget()
    };
    assert!(
        project_viscosity(&mesh, &[law.clone()], &model, &vectors, too_small, || false).is_err()
    );
    let mut calls = 0;
    let cancelled = project_viscosity(&mesh, &[law], &model, &vectors, budget(), || {
        calls += 1;
        calls >= 4
    });
    assert_eq!(cancelled, Err(fs_plate::PlateError::Cancelled));
    assert!(projection_terms(usize::MAX, 1, 2).is_none());
}
