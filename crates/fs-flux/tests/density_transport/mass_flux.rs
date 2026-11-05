//! G0/G1/G3: exact quadratic face moments, local integration by parts,
//! constant divergence and agreement with actual Euler density inventories.

use super::{Cell, Error, FluxSystem, Mesh2, TriMesh, counterflow, near};
use fs_flux::bdm::{eval_basis, tri_quad};

fn states(densities: &[f64]) -> Vec<Cell<0>> {
    densities
        .iter()
        .map(|&density_kg_m3| Cell {
            density_kg_m3,
            specific: [],
        })
        .collect()
}

fn native_value(system: &FluxSystem<'_>, velocity: &[f64], cell: usize, x: [f64; 2]) -> [f64; 2] {
    let mut value = [0.0; 2];
    for local in 0..6 {
        let edge = system.mesh.tri_edges[cell][local / 2].0;
        let basis = eval_basis(&system.bases[cell], local, x);
        for axis in 0..2 {
            value[axis] += velocity[2 * edge + local % 2] * basis[axis];
        }
    }
    value
}

#[test]
fn uniform_density_reproduces_original_full_bdm_velocity() {
    let mut geometry = Mesh2::triangles(2.0, 0.8, 2, 1);
    for p in &mut geometry.nodes {
        *p = [2.0 + p[0] + 0.4 * p[1], -3.0 + p[1]];
    }
    let mesh = TriMesh::from_mesh2(&geometry);
    let system = FluxSystem::new(&mesh);
    let mut velocity = vec![0.0; system.n_u];
    // Independent odd moments on every edge leave all cell divergences zero,
    // while exercising every P1 velocity mode and both boundary directions.
    for edge in 0..mesh.edges.len() {
        velocity[2 * edge + 1] = (edge as f64 + 1.0) * 0.03;
    }
    let rho = 3.7;
    let cells = states(&vec![rho; mesh.tris.len()]);
    let inflow: Vec<_> = mesh
        .edges
        .iter()
        .map(|edge| (edge.tris.1 == usize::MAX).then_some(cells[0]))
        .collect();
    let plan = system
        .prepare_density_transport(&velocity, 0.7, 1e-12, 1 << 20, || false)
        .unwrap();
    let flux = plan
        .prepare_mass_flux(&cells, &inflow, 1 << 20, || false)
        .unwrap();
    assert!(core::ptr::eq(flux.mesh(), &mesh));
    for cell in 0..mesh.tris.len() {
        let p = mesh.tris[cell].map(|vertex| mesh.verts[vertex]);
        for (point, _) in tri_quad(p, mesh.areas[cell]) {
            let expected = native_value(&system, &velocity, cell, point);
            let actual = flux.cell_value(cell, point).unwrap();
            for axis in 0..2 {
                near(actual[axis], rho * expected[axis], 3e-12);
            }
        }
        near(flux.cell_divergence(cell).unwrap(), 0.0, 0.0);
    }
}

fn polynomial(kind: usize, p: [f64; 2]) -> (f64, [f64; 2]) {
    match kind {
        0 => (1.0, [0.0, 0.0]),
        1 => (p[0], [1.0, 0.0]),
        2 => (p[1], [0.0, 1.0]),
        3 => (p[0] * p[0], [2.0 * p[0], 0.0]),
        4 => (p[0] * p[1], [p[1], p[0]]),
        5 => (p[1] * p[1], [0.0, 2.0 * p[1]]),
        _ => unreachable!(),
    }
}

#[test]
fn unequal_density_counterflow_closes_every_quadratic_integration_by_parts_identity() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let edge_id = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let edge = &mesh.edges[edge_id];
    let mut velocity = vec![0.0; system.n_u];
    velocity[2 * edge_id + 1] = 1.0 / 12.0; // normal velocity s-1/2
    let cells = states(&[2.0, 8.0]);
    let inflow = vec![None; mesh.edges.len()];
    let plan = system
        .prepare_density_transport(&velocity, 1.0, 1e-14, 1 << 20, || false)
        .unwrap();
    let flux = plan
        .prepare_mass_flux(&cells, &inflow, 1 << 20, || false)
        .unwrap();
    let (a, b) = (mesh.verts[edge.verts.0], mesh.verts[edge.verts.1]);
    let mut observed = [0.0; 3];
    let mut projected = [0.0; 3];
    for point in flux.edge_quadrature(edge_id).unwrap() {
        let s = ((point.point_m[0] - a[0]) * (b[0] - a[0])
            + (point.point_m[1] - a[1]) * (b[1] - a[1]))
            / edge.len.powi(2);
        for (k, p) in [1.0, 2.0 * s - 1.0, 6.0 * s * s - 6.0 * s + 1.0]
            .into_iter()
            .enumerate()
        {
            observed[k] += point.mass_rate_kg_m_s * p;
        }
    }
    // The projected trace itself is quadratic. Three-point Gauss integrates
    // its product with each quadratic test exactly, independently of the
    // original two-point-per-sign flux quadrature.
    for (s, weight) in [
        (0.5 - (3.0_f64 / 5.0).sqrt() / 2.0, 5.0 / 18.0),
        (0.5, 4.0 / 9.0),
        (0.5 + (3.0_f64 / 5.0).sqrt() / 2.0, 5.0 / 18.0),
    ] {
        let p = core::array::from_fn(|axis| a[axis] + s * (b[axis] - a[axis]));
        let left = flux.cell_value(edge.tris.0, p).unwrap();
        let right = flux.cell_value(edge.tris.1, p).unwrap();
        let normal = left[0] * edge.normal[0] + left[1] * edge.normal[1];
        near(
            right[0] * edge.normal[0] + right[1] * edge.normal[1],
            normal,
            2e-13,
        );
        for (k, polynomial) in [1.0, 2.0 * s - 1.0, 6.0 * s * s - 6.0 * s + 1.0]
            .into_iter()
            .enumerate()
        {
            projected[k] += weight * edge.len * normal * polynomial;
        }
    }
    let rho_left = cells[edge.tris.0].density_kg_m3;
    let rho_right = cells[edge.tris.1].density_kg_m3;
    let exact = [
        edge.len * (rho_left - rho_right) / 8.0,
        edge.len * (rho_left + rho_right) / 12.0,
        edge.len * (rho_left - rho_right) / 32.0,
    ];
    for k in 0..3 {
        near(observed[k], exact[k], 2e-14);
        near(projected[k], exact[k], 2e-13);
    }
    assert!(
        exact[2].abs() > 0.1,
        "BDM1 cannot retain this quadratic trace moment"
    );
    for cell in 0..mesh.tris.len() {
        let vertices = mesh.tris[cell].map(|vertex| mesh.verts[vertex]);
        for kind in 0..6 {
            let volume: f64 = tri_quad(vertices, mesh.areas[cell])
                .into_iter()
                .map(|(p, weight)| {
                    let (value, gradient) = polynomial(kind, p);
                    let j = flux.cell_value(cell, p).unwrap();
                    weight
                        * (j[0] * gradient[0]
                            + j[1] * gradient[1]
                            + flux.cell_divergence(cell).unwrap() * value)
                })
                .sum();
            let boundary: f64 = mesh.tri_edges[cell]
                .iter()
                .map(|&(edge, sign)| {
                    flux.edge_quadrature(edge)
                        .unwrap()
                        .iter()
                        .map(|q| sign * q.mass_rate_kg_m_s * polynomial(kind, q.point_m).0)
                        .sum::<f64>()
                })
                .sum();
            near(volume, boundary, 6e-13);
        }
    }
}

#[test]
fn mass_flux_divergence_matches_accepted_density_and_original_directional_transfers() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let edge = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let velocity = counterflow(&system, edge, 0.2);
    for contrast in [8.0, 1e12] {
        for depth in [0.25, 3.0] {
            let cells = states(&[2.0, contrast]);
            let inflow = vec![None; mesh.edges.len()];
            let plan = system
                .prepare_density_transport(&velocity, depth, 1e-14, 1 << 20, || false)
                .unwrap();
            let flux = plan
                .prepare_mass_flux(&cells, &inflow, 1 << 20, || false)
                .unwrap();
            let mut scratch = plan.workspace::<0>(1 << 20).unwrap();
            let mut new = cells.clone();
            let dt = 0.125;
            plan.advance_into(&cells, &inflow, dt, 0.8, &mut new, &mut scratch, || false)
                .unwrap();
            for cell in 0..cells.len() {
                near(
                    (new[cell].density_kg_m3 - cells[cell].density_kg_m3) / dt,
                    -flux.cell_divergence(cell).unwrap(),
                    2e-14 * contrast,
                );
            }
            for (edge, accepted) in scratch.face_transfers().unwrap().iter().enumerate() {
                let rates = flux.edge_directional_mass_rates(edge).unwrap();
                for direction in 0..2 {
                    near(
                        accepted.mass_kg[direction] / (dt * depth),
                        rates[direction],
                        2e-15 * contrast,
                    );
                }
            }
        }
    }
}

#[test]
fn reconstruction_is_independent_of_triangle_corner_order_and_rigid_frame() {
    let geometry = Mesh2::triangles(1.3, 0.7, 1, 1);
    let first = TriMesh::from_mesh2(&geometry);
    let first_system = FluxSystem::new(&first);
    let original_edge = first
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let first_velocity = counterflow(&first_system, original_edge, 0.2);
    let cells = states(&[2.0, 8.0]);
    let first_plan = first_system
        .prepare_density_transport(&first_velocity, 1.0, 1e-14, 1 << 20, || false)
        .unwrap();
    let first_flux = first_plan
        .prepare_mass_flux(&cells, &vec![None; first.edges.len()], 1 << 20, || false)
        .unwrap();
    let (sin, cos) = 0.43_f64.sin_cos();
    let rotate = |p: [f64; 2]| [cos * p[0] - sin * p[1], sin * p[0] + cos * p[1]];
    let mut moved = geometry;
    for p in &mut moved.nodes {
        let q = rotate(*p);
        *p = [4.0 + 3.0 * q[0], -2.0 + 3.0 * q[1]];
    }
    for tri in &mut moved.elems {
        tri.rotate_left(1);
    }
    let second = TriMesh::from_mesh2(&moved);
    let second_system = FluxSystem::new(&second);
    let mut velocity = vec![0.0; second_system.n_u];
    for (e, edge) in second.edges.iter().enumerate() {
        let old = first
            .edges
            .iter()
            .position(|old| old.verts == edge.verts)
            .unwrap();
        velocity[2 * e] = first_velocity[2 * old];
        velocity[2 * e + 1] = first_velocity[2 * old + 1];
    }
    let plan = second_system
        .prepare_density_transport(&velocity, 1.0, 1e-14, 1 << 20, || false)
        .unwrap();
    let flux = plan
        .prepare_mass_flux(&cells, &vec![None; second.edges.len()], 1 << 20, || false)
        .unwrap();
    for cell in 0..first.tris.len() {
        for (point, _) in tri_quad(first.tris[cell].map(|v| first.verts[v]), first.areas[cell]) {
            let p = rotate(point);
            let actual = flux
                .cell_value(cell, [4.0 + 3.0 * p[0], -2.0 + 3.0 * p[1]])
                .unwrap();
            let expected = rotate(first_flux.cell_value(cell, point).unwrap());
            for axis in 0..2 {
                near(actual[axis], expected[axis], 2e-12);
            }
        }
        near(
            flux.cell_divergence(cell).unwrap(),
            first_flux.cell_divergence(cell).unwrap() / 3.0,
            1e-14,
        );
    }
}

#[test]
fn preparation_refuses_missing_inflow_invalid_state_storage_and_late_cancellation() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let mut velocity = vec![0.0; system.n_u];
    for e in 0..mesh.edges.len() {
        velocity[2 * e] = mesh.edges[e].normal[0];
    }
    let plan = system
        .prepare_density_transport(&velocity, 1.0, 1e-12, 1 << 20, || false)
        .unwrap();
    let cells = states(&[2.0, 8.0]);
    let missing = vec![None; mesh.edges.len()];
    assert!(matches!(
        plan.prepare_mass_flux(&cells, &missing, 1 << 20, || false),
        Err(Error::Invalid(_))
    ));
    let inflow: Vec<_> = mesh
        .edges
        .iter()
        .map(|edge| {
            (edge.tris.1 == usize::MAX).then_some(Cell {
                density_kg_m3: 4.0,
                specific: [],
            })
        })
        .collect();
    assert!(matches!(
        plan.prepare_mass_flux(&cells, &inflow, 0, || false),
        Err(Error::Storage { .. })
    ));
    let mut invalid = cells.clone();
    invalid[0].density_kg_m3 = -1.0;
    assert!(matches!(
        plan.prepare_mass_flux(&invalid, &inflow, 1 << 20, || false),
        Err(Error::Invalid(_))
    ));
    let mut calls = 0;
    assert!(matches!(
        plan.prepare_mass_flux(&cells, &inflow, 1 << 20, || {
            calls += 1;
            calls > 1 + cells.len() + mesh.edges.len()
        }),
        Err(Error::Cancelled)
    ));
    let flux = plan
        .prepare_mass_flux(&cells, &inflow, 1 << 20, || false)
        .unwrap();
    assert!(flux.cell_value(usize::MAX, [0.0; 2]).is_err());
    assert!(flux.cell_value(0, [f64::NAN, 0.0]).is_err());
    assert!(flux.edge_quadrature(usize::MAX).is_none());
    assert_eq!(cells, states(&[2.0, 8.0]));
}
