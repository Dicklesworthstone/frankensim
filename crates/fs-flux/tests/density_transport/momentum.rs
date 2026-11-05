//! Conservative native momentum, changing density and kinetic-energy balance.

use super::{Cell, Error, FluxSystem, Mesh2, TriMesh, near};
use fs_flux::FluxParams;
use fs_flux::bdm::{edge_gauss_pub, eval_basis, tri_quad};

fn states(mesh: &TriMesh, uniform: Option<f64>) -> Vec<Cell<0>> {
    (0..mesh.tris.len())
        .map(|cell| Cell {
            density_kg_m3: uniform.unwrap_or(2.0 + 3.0 * cell as f64),
            specific: [],
        })
        .collect()
}

fn advector(system: &FluxSystem<'_>, exterior: bool) -> Vec<f64> {
    let mut field = vec![0.0; system.n_u];
    for (edge, geometry) in system.mesh.edges.iter().enumerate() {
        if exterior || geometry.tris.1 != usize::MAX {
            // Zero mean on every edge guarantees cellwise volume continuity;
            // independent linear moments produce true reversing face streams.
            field[2 * edge + 1] = 0.012 * (edge + 1) as f64;
        }
    }
    field
}

fn value(system: &FluxSystem<'_>, field: &[f64], cell: usize, point: [f64; 2]) -> [f64; 2] {
    let mut value = [0.0; 2];
    for local in 0..6 {
        let coefficient = field[2 * system.mesh.tri_edges[cell][local / 2].0 + local % 2];
        let basis = eval_basis(&system.bases[cell], local, point);
        for axis in 0..2 {
            value[axis] += coefficient * basis[axis];
        }
    }
    value
}

fn exterior(point: [f64; 2]) -> [f64; 2] {
    [0.3 + 0.7 * point[0] - 0.2 * point[1], -0.8 + 0.4 * point[1]]
}

fn samples(mesh: &TriMesh) -> Vec<Option<[[f64; 2]; 2]>> {
    mesh.edges
        .iter()
        .map(|edge| {
            (edge.tris.1 == usize::MAX).then(|| {
                edge_gauss_pub(mesh.verts[edge.verts.0], mesh.verts[edge.verts.1])
                    .map(|(point, _)| exterior(point))
            })
        })
        .collect()
}

fn constant(mesh: &TriMesh, vector: [f64; 2]) -> Vec<f64> {
    let mut field = vec![0.0; 2 * mesh.edges.len()];
    for (index, edge) in mesh.edges.iter().enumerate() {
        field[2 * index] = vector[0] * edge.normal[0] + vector[1] * edge.normal[1];
    }
    field
}

fn norm2(value: [f64; 2]) -> f64 {
    value[0] * value[0] + value[1] * value[1]
}

#[test]
fn uniform_density_matches_original_native_convection_including_reversing_inflow() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.3, 0.8, 2, 2));
    let system = FluxSystem::new(&mesh);
    let rho = 3.7;
    let density = states(&mesh, Some(rho));
    let density_inflow: Vec<_> = mesh
        .edges
        .iter()
        .map(|edge| (edge.tris.1 == usize::MAX).then_some(density[0]))
        .collect();
    let velocity: Vec<_> = (0..system.n_u).map(|i| (0.3 * i as f64).sin()).collect();
    let w = advector(&system, true);
    let plan = system
        .prepare_density_transport(&w, 1.0, 1e-12, 1 << 20, || false)
        .unwrap();
    let flux = plan
        .prepare_mass_flux(&density, &density_inflow, 1 << 20, || false)
        .unwrap();
    let mut scratch = flux.momentum_workspace(1 << 20).unwrap();
    let mut action = vec![0.0; system.n_u];
    flux.apply_momentum_into(
        &system,
        &velocity,
        &samples(&mesh),
        &mut action,
        &mut scratch,
        || false,
    )
    .unwrap();

    // Original production assembly is an independent volume-velocity oracle.
    // Subtraction removes pressure coupling and prescribed-normal rows.
    let params = FluxParams {
        nu: 0.0,
        ..FluxParams::default()
    };
    let (transport, load) = system.assemble(params, &|_| [0.0; 2], &exterior, Some(&w));
    let (base, base_load) = system.assemble(params, &|_| [0.0; 2], &exterior, None);
    let mut full = velocity.clone();
    full.resize(system.n, 0.0);
    let mut left = vec![0.0; system.n];
    let mut right = vec![0.0; system.n];
    transport.spmv(&full, &mut left);
    base.spmv(&full, &mut right);
    let mut checked = 0;
    for (edge, geometry) in mesh.edges.iter().enumerate() {
        if geometry.tris.1 != usize::MAX {
            for row in [2 * edge, 2 * edge + 1] {
                let expected = rho * (left[row] - right[row] - load[row] + base_load[row]);
                near(action[row], expected, 3e-11 * expected.abs().max(1.0));
                checked += 1;
            }
        }
    }
    assert!(checked > 4);
}

#[test]
fn transported_density_and_native_momentum_preserve_a_constant_velocity() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 0.7, 2, 1));
    let system = FluxSystem::new(&mesh);
    let density = states(&mesh, None);
    let boundary = vec![None; mesh.edges.len()];
    let w = advector(&system, false);
    let plan = system
        .prepare_density_transport(&w, 0.4, 1e-12, 1 << 20, || false)
        .unwrap();
    let flux = plan
        .prepare_mass_flux(&density, &boundary, 1 << 20, || false)
        .unwrap();
    let dt = 0.02;
    let mut advanced = density.clone();
    let mut density_scratch = plan.workspace::<0>(1 << 20).unwrap();
    plan.advance_into(
        &density,
        &boundary,
        dt,
        0.8,
        &mut advanced,
        &mut density_scratch,
        || false,
    )
    .unwrap();
    assert!(advanced.iter().zip(&density).any(|(new, old)| new != old));
    let uniform = [1.3, -0.8];
    let velocity = constant(&mesh, uniform);
    let mut action = vec![0.0; system.n_u];
    let mut scratch = flux.momentum_workspace(1 << 20).unwrap();
    flux.apply_momentum_into(
        &system,
        &velocity,
        &vec![None; mesh.edges.len()],
        &mut action,
        &mut scratch,
        || false,
    )
    .unwrap();

    // Assemble the ACTUAL changing momentum inventory, not a projected
    // velocity or a hand-specified density derivative. Every native test row
    // must cancel convection when old/new velocity is the same constant.
    let mut inertia = vec![0.0; system.n_u];
    for cell in 0..mesh.tris.len() {
        let rho_dot = (advanced[cell].density_kg_m3 - density[cell].density_kg_m3) / dt;
        let p = mesh.tris[cell].map(|vertex| mesh.verts[vertex]);
        for (point, weight) in tri_quad(p, mesh.areas[cell]) {
            for local in 0..6 {
                let test = eval_basis(&system.bases[cell], local, point);
                let row = 2 * mesh.tri_edges[cell][local / 2].0 + local % 2;
                inertia[row] += weight * rho_dot * (uniform[0] * test[0] + uniform[1] * test[1]);
            }
        }
    }
    for row in 0..system.n_u {
        near(inertia[row] + action[row], 0.0, 2e-11);
    }
}

#[test]
fn unequal_density_native_momentum_closes_energy_and_boundary_momentum() {
    let mut geometry = Mesh2::triangles(1.1, 0.6, 2, 2);
    for p in &mut geometry.nodes {
        *p = [2.0 + p[0] + 0.3 * p[1], -1.0 + p[1]];
    }
    let mesh = TriMesh::from_mesh2(&geometry);
    let system = FluxSystem::new(&mesh);
    let density = states(&mesh, None);
    let exterior_density: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .map(|(index, edge)| {
            (edge.tris.1 == usize::MAX).then_some(Cell {
                density_kg_m3: 4.0 + index as f64,
                specific: [],
            })
        })
        .collect();
    let w = advector(&system, true);
    let plan = system
        .prepare_density_transport(&w, 1.0, 1e-12, 1 << 20, || false)
        .unwrap();
    let flux = plan
        .prepare_mass_flux(&density, &exterior_density, 1 << 20, || false)
        .unwrap();
    let velocity: Vec<_> = (0..system.n_u).map(|i| (0.7 * i as f64).cos()).collect();
    let mut scratch = flux.momentum_workspace(1 << 20).unwrap();
    let mut action = vec![0.0; system.n_u];
    flux.apply_momentum_into(
        &system,
        &velocity,
        &samples(&mesh),
        &mut action,
        &mut scratch,
        || false,
    )
    .unwrap();
    let power: f64 = action.iter().zip(&velocity).map(|(r, u)| r * u).sum();
    let mut density_power = 0.0;
    for cell in 0..mesh.tris.len() {
        let p = mesh.tris[cell].map(|vertex| mesh.verts[vertex]);
        for (point, weight) in tri_quad(p, mesh.areas[cell]) {
            density_power += 0.5
                * weight
                * flux.cell_divergence(cell).unwrap()
                * norm2(value(&system, &velocity, cell, point));
        }
    }
    let mut boundary_power = 0.0;
    let mut loss = 0.0;
    let mut boundary_momentum = [0.0; 2];
    for (edge_id, edge) in mesh.edges.iter().enumerate() {
        for point in flux.edge_quadrature(edge_id).unwrap() {
            let rate = point.mass_rate_kg_m_s;
            let first = value(&system, &velocity, edge.tris.0, point.point_m);
            let second = if edge.tris.1 == usize::MAX {
                exterior(point.point_m)
            } else {
                value(&system, &velocity, edge.tris.1, point.point_m)
            };
            if edge.tris.1 == usize::MAX {
                let upwind = if rate >= 0.0 { first } else { second };
                boundary_power += 0.5 * rate * norm2(upwind);
                for axis in 0..2 {
                    boundary_momentum[axis] += rate * upwind[axis];
                }
            }
            if edge.tris.1 != usize::MAX || rate < 0.0 {
                loss += 0.5 * rate.abs() * norm2([first[0] - second[0], first[1] - second[1]]);
            }
        }
    }
    assert!(loss > 1.0 && density_power.abs() > 1e-3);
    near(
        power,
        density_power + boundary_power + loss,
        5e-11 * power.abs().max(1.0),
    );
    for axis in 0..2 {
        let mut unit = [0.0; 2];
        unit[axis] = 1.0;
        let test = constant(&mesh, unit);
        let total: f64 = action.iter().zip(test).map(|(r, v)| r * v).sum();
        near(total, boundary_momentum[axis], 3e-11 * total.abs().max(1.0));
    }
}

#[test]
fn native_momentum_refusals_and_final_cancellation_preserve_output_and_retry() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let cells = states(&mesh, None);
    let inlet: Vec<_> = mesh
        .edges
        .iter()
        .map(|edge| (edge.tris.1 == usize::MAX).then_some(cells[0]))
        .collect();
    let w = advector(&system, true);
    let plan = system
        .prepare_density_transport(&w, 1.0, 1e-12, 1 << 20, || false)
        .unwrap();
    let flux = plan
        .prepare_mass_flux(&cells, &inlet, 1 << 20, || false)
        .unwrap();
    assert!(matches!(
        flux.momentum_workspace(0),
        Err(Error::Storage { .. })
    ));
    let mut scratch = flux.momentum_workspace(1 << 20).unwrap();
    let velocity = constant(&mesh, [0.7, -1.1]);
    let incoming = samples(&mesh);
    let mut output = vec![123.0; system.n_u];
    assert!(matches!(
        flux.apply_momentum_into(
            &system,
            &velocity,
            &vec![None; mesh.edges.len()],
            &mut output,
            &mut scratch,
            || false
        ),
        Err(Error::Invalid(_))
    ));
    assert_eq!(output, vec![123.0; system.n_u]);
    let other_mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let other_system = FluxSystem::new(&other_mesh);
    assert!(matches!(
        flux.apply_momentum_into(
            &other_system,
            &velocity,
            &incoming,
            &mut output,
            &mut scratch,
            || false
        ),
        Err(Error::Invalid(_))
    ));
    assert_eq!(output, vec![123.0; system.n_u]);
    let mut polls = 0;
    flux.apply_momentum_into(
        &system,
        &velocity,
        &incoming,
        &mut output,
        &mut scratch,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let expected = output.clone();
    output.fill(456.0);
    let mut seen = 0;
    assert_eq!(
        flux.apply_momentum_into(
            &system,
            &velocity,
            &incoming,
            &mut output,
            &mut scratch,
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(Error::Cancelled)
    );
    assert_eq!(seen, polls);
    assert_eq!(output, vec![456.0; system.n_u]);
    flux.apply_momentum_into(
        &system,
        &velocity,
        &incoming,
        &mut output,
        &mut scratch,
        || false,
    )
    .unwrap();
    assert_eq!(output, expected);
}
