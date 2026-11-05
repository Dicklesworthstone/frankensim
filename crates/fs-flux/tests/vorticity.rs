//! Native BDM curl, distributional interface sheets and solved-flow accuracy.
#![cfg(feature = "continuum")]

use fs_flux::bdm::{edge_gauss_pub, eval_basis};
use fs_flux::ns::kinematics::FluxKinematicsError;
use fs_flux::{FluxParams, FluxSystem, TriMesh};
use fs_solid::Mesh2;

fn interpolate(system: &FluxSystem<'_>, field: impl Fn(usize, [f64; 2]) -> [f64; 2]) -> Vec<f64> {
    let mesh = system.mesh;
    let mut coefficients = vec![0.0; system.n_u];
    for (index, edge) in mesh.edges.iter().enumerate() {
        let (a, b) = (mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]);
        for (point, weight) in edge_gauss_pub(a, b) {
            let value = field(edge.tris.0, point);
            let normal = value[0] * edge.normal[0] + value[1] * edge.normal[1];
            let coordinate = ((point[0] - a[0]) * (b[0] - a[0])
                + (point[1] - a[1]) * (b[1] - a[1]))
                / edge.len.powi(2)
                - 0.5;
            coefficients[2 * index] += weight * normal / edge.len;
            coefficients[2 * index + 1] += weight * normal * coordinate / edge.len;
        }
    }
    coefficients
}

fn local_velocity(system: &FluxSystem<'_>, field: &[f64], cell: usize, p: [f64; 2]) -> [f64; 2] {
    let mut velocity = [0.0; 2];
    for local in 0..6 {
        let (edge, _) = system.mesh.tri_edges[cell][local / 2];
        let basis = eval_basis(&system.bases[cell], local, p);
        for component in 0..2 {
            velocity[component] += field[2 * edge + local % 2] * basis[component];
        }
    }
    velocity
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn affine_kinematics_distinguish_rotation_shear_and_irrotational_strain() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(3.0, 0.7, 3, 2));
    let system = FluxSystem::new(&mesh);
    for (gradient, curl, divergence, shear) in [
        ([[0.0, -1.7], [1.7, 0.0]], 3.4, 0.0, 0.0),
        ([[0.0, 2.5], [0.0, 0.0]], -2.5, 0.0, 2.5),
        ([[1.2, 0.0], [0.0, -1.2]], 0.0, 0.0, 2.4),
        ([[0.4, 0.0], [0.0, 0.6]], 0.0, 1.0, 1.04_f64.sqrt()),
    ] {
        let mut field = interpolate(&system, |_, p| {
            [
                0.31 + gradient[0][0] * p[0] + gradient[0][1] * p[1],
                -0.21 + gradient[1][0] * p[0] + gradient[1][1] * p[1],
            ]
        });
        // Pressure entries play no role in physical velocity diagnostics.
        field.extend(vec![f64::NAN; mesh.tris.len()]);
        for cell in system.cell_kinematics(&field).unwrap() {
            for (actual, expected) in cell
                .gradient()
                .iter()
                .flatten()
                .zip(gradient.iter().flatten())
            {
                close(*actual, *expected, 3e-12);
            }
            close(cell.vorticity(), curl, 3e-12);
            close(cell.divergence(), divergence, 3e-12);
            close(cell.shear_rate(), shear, 3e-12);
        }
        let balance = system.vorticity_balance(&field).unwrap();
        close(balance.cell_circulation, curl * 2.1, 3e-11);
        close(balance.sheet_circulation, 0.0, 3e-11);
        close(balance.boundary_circulation, curl * 2.1, 3e-11);
        close(balance.stokes_defect(), 0.0, 3e-11);
    }
}

fn two_sided_mesh() -> TriMesh {
    TriMesh::from_mesh2(&Mesh2 {
        nodes: vec![[0.0, 0.0], [2.0, 0.0], [1.0, 1.0], [1.0, -1.0]],
        elems: vec![vec![0, 1, 2], vec![1, 0, 3]],
        patches: vec![],
    })
}

#[test]
fn uniform_tangential_jump_has_vorticity_even_when_every_cell_curl_is_zero() {
    let mesh = two_sided_mesh();
    let system = FluxSystem::new(&mesh);
    let field = interpolate(&system, |cell, _| [if cell == 0 { 0.5 } else { 2.0 }, 0.0]);
    let edge = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let sheet = system.vorticity_sheet(&field, edge).unwrap().unwrap();
    close(sheet.endpoint_strengths[0], 1.5, 1e-13);
    close(sheet.endpoint_strengths[1], 1.5, 1e-13);
    close(sheet.circulation(), 3.0, 1e-13);
    let balance = system.vorticity_balance(&field).unwrap();
    close(balance.cell_circulation, 0.0, 1e-13);
    close(balance.sheet_circulation, 3.0, 1e-13);
    close(balance.boundary_circulation, 3.0, 1e-13);
    close(balance.stokes_defect(), 0.0, 1e-13);
    close(
        system.weak_vorticity(&field, &[1.0; 4]).unwrap(),
        3.0,
        1e-13,
    );
}

#[test]
fn affine_sheet_with_zero_circulation_retains_its_weak_first_moment() {
    let mesh = two_sided_mesh();
    let system = FluxSystem::new(&mesh);
    // Both traces have zero normal velocity along y=0; each cell is div-free.
    let field = interpolate(&system, |cell, p| {
        if cell == 0 {
            [0.0; 2]
        } else {
            [p[0] - 1.0, -p[1]]
        }
    });
    let edge = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let sheet = system.vorticity_sheet(&field, edge).unwrap().unwrap();
    close(sheet.endpoint_strengths[0], -1.0, 1e-13);
    close(sheet.endpoint_strengths[1], 1.0, 1e-13);
    close(sheet.circulation(), 0.0, 1e-13);
    for cell in system.cell_kinematics(&field).unwrap() {
        close(cell.vorticity(), 0.0, 1e-13);
        close(cell.divergence(), 0.0, 1e-13);
    }
    let scalar = mesh.verts.iter().map(|p| p[0]).collect::<Vec<_>>();
    // Integral_0^2 (x - 1) x dx = 2/3, despite zero total circulation.
    close(
        system.weak_vorticity(&field, &scalar).unwrap(),
        2.0 / 3.0,
        1e-13,
    );
}

#[test]
fn arbitrary_bdm_field_obeys_stokes_and_nonconstant_weak_integration_by_parts() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.7, 0.9, 4, 3));
    let system = FluxSystem::new(&mesh);
    let field = (0..system.n_u)
        .map(|i| (i as f64 * 0.73 + 0.11).sin())
        .collect::<Vec<_>>();
    let balance = system.vorticity_balance(&field).unwrap();
    assert!(balance.sheet_circulation.abs() > 0.1);
    close(balance.stokes_defect(), 0.0, 4e-12);

    // Independent RHS: boundary phi*u.t + integral(u_x*phi_y-u_y*phi_x).
    let scalar = |p: [f64; 2]| 0.7 + 0.4 * p[0] - 1.3 * p[1];
    let mut expected = 0.0;
    for (cell, &centroid) in mesh.centroids.iter().enumerate() {
        let u = local_velocity(&system, &field, cell, centroid);
        expected += mesh.areas[cell] * (-1.3 * u[0] - 0.4 * u[1]);
    }
    for edge in &mesh.edges {
        if edge.tris.1 != usize::MAX {
            continue;
        }
        for (point, weight) in edge_gauss_pub(mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]) {
            let u = local_velocity(&system, &field, edge.tris.0, point);
            expected += weight * scalar(point) * (u[1] * edge.normal[0] - u[0] * edge.normal[1]);
        }
    }
    let values = mesh.verts.iter().copied().map(scalar).collect::<Vec<_>>();
    close(
        system.weak_vorticity(&field, &values).unwrap(),
        expected,
        4e-12,
    );
}

#[test]
fn solved_poiseuille_vorticity_converges_under_mesh_refinement() {
    let exact = |p: [f64; 2]| [p[1] * (1.0 - p[1]), 0.0];
    let mut errors = Vec::new();
    for subdivisions in [2, 4, 8] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, subdivisions, subdivisions));
        let system = FluxSystem::new(&mesh);
        let parameters = FluxParams::default();
        let (operator, rhs) = system.assemble(parameters, &|_| [2.0, 0.0], &exact, None);
        let solution = system.solve_linear(&operator, &rhs, parameters);
        assert!(solution.converged(parameters.tol));
        let cells = system.cell_kinematics(&solution.x).unwrap();
        let mut error2 = 0.0;
        for (cell, kinematics) in cells.iter().enumerate() {
            close(kinematics.divergence(), 0.0, 1e-9);
            // Exact integral of (constant native curl - (2y-1))² on a triangle.
            let residuals =
                mesh.tris[cell].map(|v| kinematics.vorticity() - (2.0 * mesh.verts[v][1] - 1.0));
            error2 += mesh.areas[cell] / 6.0
                * (residuals.iter().map(|r| r * r).sum::<f64>()
                    + residuals[0] * residuals[1]
                    + residuals[0] * residuals[2]
                    + residuals[1] * residuals[2]);
        }
        let balance = system.vorticity_balance(&solution.x).unwrap();
        close(balance.stokes_defect(), 0.0, 1e-11);
        errors.push(error2.sqrt());
    }
    println!("solved Poiseuille native-curl L2 errors (2, 4, 8 subdivisions): {errors:?}");
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 1.7, "native curl L2 errors: {errors:?}");
    }
}

#[test]
fn checked_queries_reject_incomplete_nonfinite_and_out_of_range_data() {
    let mesh = two_sided_mesh();
    let system = FluxSystem::new(&mesh);
    let mut field = vec![0.0; system.n_u];
    assert!(system.cell_kinematics(&field[..field.len() - 1]).is_err());
    assert!(system.try_cell_kinematics(&field, 2).is_none());
    assert!(system.vorticity_sheet(&field, mesh.edges.len()).is_err());
    assert!(system.weak_vorticity(&field, &[1.0; 3]).is_err());
    assert!(system.weak_vorticity(&field, &[f64::NAN; 4]).is_err());
    let boundary = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 == usize::MAX)
        .unwrap();
    assert_eq!(system.vorticity_sheet(&field, boundary).unwrap(), None);
    field[0] = f64::NAN;
    assert!(matches!(
        system.cell_kinematics(&field),
        Err(FluxKinematicsError::NonFiniteCell(_))
    ));
    assert!(system.vorticity_balance(&field).is_err());
}
