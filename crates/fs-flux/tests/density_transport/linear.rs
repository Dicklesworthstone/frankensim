//! Bounded spatial density reconstruction and its actual accepted mass flux.

use super::{Cell, Error, FluxSystem, Mesh2, TriMesh, counterflow, near};
use fs_flux::bdm::tri_quad;
use fs_flux::ns::transport::density::DensityReconstruction;
use fs_flux::ns::transport::density::ssprk2::DensitySspRk2Controls;
use fs_flux::ns::transport::density::subcycle::DensitySubcycleControls;

fn uniform_velocity(system: &FluxSystem<'_>, velocity: [f64; 2]) -> Vec<f64> {
    let mut moments = vec![0.0; system.n_u];
    for (edge, face) in system.mesh.edges.iter().enumerate() {
        moments[2 * edge] = velocity[0] * face.normal[0] + velocity[1] * face.normal[1];
    }
    moments
}

fn profile(x: f64, time: f64) -> f64 {
    1.0 + 0.8 * (-((x - 0.4 - 0.5 * time) / 0.09).powi(2)).exp()
}

fn profile_average(mesh: &TriMesh, cell: usize, time: f64) -> f64 {
    tri_quad(
        mesh.tris[cell].map(|vertex| mesh.verts[vertex]),
        mesh.areas[cell],
    )
    .into_iter()
    .map(|(point, weight)| weight * profile(point[0], time))
    .sum::<f64>()
        / mesh.areas[cell]
}

fn translation_error(nx: usize, reconstruction: DensityReconstruction) -> f64 {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, nx, 4));
    let system = FluxSystem::new(&mesh).with_density_reconstruction(reconstruction);
    let plan = system
        .prepare_density_transport(
            &uniform_velocity(&system, [0.5, 0.0]),
            1.0,
            1e-12,
            1 << 22,
            || false,
        )
        .unwrap();
    let initial: Vec<_> = (0..mesh.tris.len())
        .map(|cell| Cell {
            density_kg_m3: profile_average(&mesh, cell, 0.0),
            specific: [123.0],
        })
        .collect();
    // The pulse is far from both open ends: its omitted inlet perturbation
    // is below 3e-9. The comparison uses actual triangle-average analytic data.
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|edge| {
            mesh.is_boundary(edge).then_some(Cell {
                density_kg_m3: 1.0,
                specific: [123.0],
            })
        })
        .collect();
    let mut output = initial.clone();
    let mut work = plan.ssprk2_workspace::<1>(1 << 22).unwrap();
    // Match BOTH schemes' physical time subdivisions to isolate spatial error.
    let maximum_courant = if reconstruction == DensityReconstruction::BoundedLinear {
        0.6
    } else {
        0.2
    };
    let report = work
        .advance_into(
            &initial,
            &inflow,
            0.3,
            DensitySspRk2Controls {
                maximum_courant,
                maximum_substeps: 256,
            },
            &mut output,
            || false,
        )
        .unwrap();
    assert!(report.maximum_stage_courant <= 0.2);
    near(report.transport.relative_mass_balance_residual, 0.0, 2e-13);
    assert!(output.iter().all(|cell| cell.specific == [123.0]));
    output
        .iter()
        .enumerate()
        .map(|(cell, state)| {
            mesh.areas[cell] * (state.density_kg_m3 - profile_average(&mesh, cell, 0.3)).abs()
        })
        .sum()
}

#[test]
fn bounded_linear_density_reduces_smooth_translation_error_with_refinement() {
    let donor = translation_error(64, DensityReconstruction::DonorCell);
    let coarse = translation_error(32, DensityReconstruction::BoundedLinear);
    let fine = translation_error(64, DensityReconstruction::BoundedLinear);
    assert!(fine < 0.6 * donor, "linear={fine:e}, donor={donor:e}");
    assert!(coarse / fine > 2.1, "coarse={coarse:e}, fine={fine:e}");
}

#[test]
fn affine_density_counterflow_uses_each_streams_exact_flux_centroid() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 4, 4));
    let system =
        FluxSystem::new(&mesh).with_density_reconstruction(DensityReconstruction::BoundedLinear);
    let edge = mesh
        .edges
        .iter()
        .position(|face| {
            face.tris.1 != usize::MAX
                && [face.tris.0, face.tris.1].into_iter().all(|cell| {
                    mesh.tris[cell]
                        .into_iter()
                        .all(|vertex| mesh.verts[vertex].into_iter().all(|x| x > 0.0 && x < 1.0))
                })
        })
        .unwrap();
    let plan = system
        .prepare_density_transport(
            &counterflow(&system, edge, 0.02),
            1.7,
            1e-12,
            1 << 22,
            || false,
        )
        .unwrap();
    let density = |point: [f64; 2]| 2.0 + point[0] + 2.0 * point[1];
    let cells: Vec<Cell<0>> = mesh
        .centroids
        .iter()
        .map(|&point| Cell {
            density_kg_m3: density(point),
            specific: [],
        })
        .collect();
    let inflow = vec![None; mesh.edges.len()];
    let mut output = cells.clone();
    let mut work = plan.workspace::<0>(1 << 22).unwrap();
    let dt = 0.8 * plan.maximum_timestep_s(0.6).unwrap();
    plan.advance_into(&cells, &inflow, dt, 0.6, &mut output, &mut work, || false)
        .unwrap();
    let face = &mesh.edges[edge];
    let a = mesh.verts[face.verts.0];
    let b = mesh.verts[face.verts.1];
    let rates = plan.edge_volume_rates()[edge];
    let transfer = work.face_transfers().unwrap()[edge];
    let flux = plan
        .prepare_mass_flux(&cells, &inflow, 1 << 22, || false)
        .unwrap();
    for (direction, s, rate) in [
        (0, 5.0 / 6.0, rates.along_normal_m2_s()),
        (1, 1.0 / 6.0, rates.against_normal_m2_s()),
    ] {
        let point = core::array::from_fn(|axis| a[axis] + s * (b[axis] - a[axis]));
        let expected = dt * 1.7 * rate * density(point);
        near(transfer.mass_kg[direction], expected, 1e-13);
        near(
            dt * 1.7 * flux.edge_directional_mass_rates(edge).unwrap()[direction],
            expected,
            1e-13,
        );
    }
    for cell in 0..cells.len() {
        near(
            output[cell].density_kg_m3 - cells[cell].density_kg_m3
                + dt * flux.cell_divergence(cell).unwrap(),
            0.0,
            2e-12,
        );
    }
    let accepted = output.clone();
    assert!(matches!(
        plan.advance_into(
            &cells,
            &inflow,
            1.1 * plan.maximum_timestep_s(0.6).unwrap(),
            0.6,
            &mut output,
            &mut work,
            || false,
        ),
        Err(Error::Courant { .. })
    ));
    assert_eq!(output, accepted);
    assert!(work.face_transfers().is_none());
}

#[test]
fn density_contrast_stays_positive_and_stage_flux_matches_every_accepted_transfer() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 8, 3));
    let system =
        FluxSystem::new(&mesh).with_density_reconstruction(DensityReconstruction::BoundedLinear);
    let depth = 2.5;
    let dt = 0.25;
    let plan = system
        .prepare_density_transport(
            &uniform_velocity(&system, [0.5, 0.0]),
            depth,
            1e-12,
            1 << 22,
            || false,
        )
        .unwrap();
    let cells: Vec<_> = mesh
        .centroids
        .iter()
        .map(|point| Cell {
            density_kg_m3: if point[0] < 0.5 { 1000.0 } else { 1.2 },
            specific: [if point[0] < 0.5 { 0.7 } else { 0.2 }, 123.0],
        })
        .collect();
    let inflow: Vec<_> = (0..mesh.edges.len())
        .map(|edge| {
            mesh.is_boundary(edge).then_some(Cell {
                density_kg_m3: 1000.0,
                specific: [0.7, 123.0],
            })
        })
        .collect();
    for rk in [false, true] {
        let mut output = cells.clone();
        let (report, transfers, flux) = if rk {
            let mut work = plan.ssprk2_workspace::<2>(1 << 22).unwrap();
            let report = work
                .advance_into(
                    &cells,
                    &inflow,
                    dt,
                    DensitySspRk2Controls {
                        maximum_courant: 0.75,
                        maximum_substeps: 128,
                    },
                    &mut output,
                    || false,
                )
                .unwrap();
            (
                report.transport,
                work.face_transfers().unwrap().to_vec(),
                work.prepare_mass_flux(1 << 22, || false).unwrap(),
            )
        } else {
            let mut work = plan.subcycle_workspace::<2>(1 << 22).unwrap();
            let report = work
                .advance_into(
                    &cells,
                    &inflow,
                    dt,
                    DensitySubcycleControls {
                        maximum_courant: 0.75,
                        maximum_substeps: 128,
                    },
                    &mut output,
                    || false,
                )
                .unwrap();
            (
                report.transport,
                work.face_transfers().unwrap().to_vec(),
                work.prepare_mass_flux(1 << 22, || false).unwrap(),
            )
        };
        for cell in &output {
            assert!((1.2 - 1e-10..=1000.0 + 1e-10).contains(&cell.density_kg_m3));
            assert!((0.2 - 1e-13..=0.7 + 1e-13).contains(&cell.specific[0]));
            assert_eq!(cell.specific[1], 123.0);
        }
        near(report.relative_mass_balance_residual, 0.0, 3e-14);
        for residual in report.relative_carried_balance_residual {
            near(residual, 0.0, 5e-14);
        }
        for (edge, transfer) in transfers.iter().enumerate() {
            let rates = flux.edge_directional_mass_rates(edge).unwrap();
            for direction in 0..2 {
                near(
                    dt * depth * rates[direction],
                    transfer.mass_kg[direction],
                    3e-13 * transfer.mass_kg[direction].max(1.0),
                );
            }
        }
        for cell in 0..cells.len() {
            near(
                output[cell].density_kg_m3 - cells[cell].density_kg_m3
                    + dt * flux.cell_divergence(cell).unwrap(),
                0.0,
                2e-10,
            );
        }
    }
}

fn thermal_average(mesh: &TriMesh, cell: usize, time: f64) -> Cell<3> {
    let mut mass = 0.0;
    let mut enthalpy = 0.0;
    for (point, weight) in tri_quad(
        mesh.tris[cell].map(|vertex| mesh.verts[vertex]),
        mesh.areas[cell],
    ) {
        let rho = profile(point[0], time);
        let h = 200_000.0 + 100_000.0 * (profile(point[0] - 0.025, time) - 1.0);
        mass += weight * rho;
        enthalpy += weight * rho * h;
    }
    Cell {
        density_kg_m3: mass / mesh.areas[cell],
        specific: [enthalpy / mass, 0.3, -123.0],
    }
}

fn thermal_translation_error(nx: usize, reconstruction: DensityReconstruction) -> f64 {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, nx, 4));
    let system = FluxSystem::new(&mesh).with_density_reconstruction(reconstruction);
    let plan = system
        .prepare_density_transport(
            &uniform_velocity(&system, [0.5, 0.0]),
            1.0,
            1e-12,
            1 << 22,
            || false,
        )
        .unwrap();
    let cells: Vec<_> = (0..mesh.tris.len())
        .map(|cell| thermal_average(&mesh, cell, 0.0))
        .collect();
    let inflow: Vec<_> = mesh
        .edges
        .iter()
        .map(|edge| {
            (edge.tris.1 == usize::MAX).then_some(Cell {
                density_kg_m3: 1.0,
                specific: [200_000.0, 0.3, -123.0],
            })
        })
        .collect();
    let mut output = cells.clone();
    let mut work = plan.ssprk2_workspace::<3>(1 << 22).unwrap();
    let report = work
        .advance_into(
            &cells,
            &inflow,
            0.3,
            DensitySspRk2Controls {
                maximum_courant: 0.6,
                maximum_substeps: 256,
            },
            &mut output,
            || false,
        )
        .unwrap();
    near(report.transport.relative_mass_balance_residual, 0.0, 3e-13);
    for residual in report.transport.relative_carried_balance_residual {
        near(residual, 0.0, 3e-13);
    }
    assert!(
        output
            .iter()
            .all(|state| state.specific[1..] == [0.3, -123.0])
    );
    output
        .iter()
        .enumerate()
        .map(|(cell, state)| {
            mesh.areas[cell]
                * (state.specific[0] - thermal_average(&mesh, cell, 0.3).specific[0]).abs()
        })
        .sum::<f64>()
        / 100_000.0
}

#[test]
fn bounded_inventory_convection_reduces_thermal_diffusion_and_refines_in_space() {
    // Identical density reconstruction, fixed velocity, physical timestep and
    // exact mass-weighted triangle-average reference isolate the carried heat.
    let donor = thermal_translation_error(64, DensityReconstruction::BoundedLinear);
    let coarse = thermal_translation_error(32, DensityReconstruction::BoundedLinearInventories);
    let fine = thermal_translation_error(64, DensityReconstruction::BoundedLinearInventories);
    assert!(
        fine < 0.55 * donor,
        "reconstructed={fine:e}, donor={donor:e}"
    );
    assert!(coarse / fine > 3.0, "coarse={coarse:e}, fine={fine:e}");
}

#[test]
fn reconstructed_counterflow_reproduces_affine_mass_means_and_cubic_inventory_flux() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 8, 8));
    let system = FluxSystem::new(&mesh)
        .with_density_reconstruction(DensityReconstruction::BoundedLinearInventories);
    let edge = mesh
        .edges
        .iter()
        .position(|face| {
            face.tris.1 != usize::MAX
                && [face.tris.0, face.tris.1].into_iter().all(|cell| {
                    mesh.tris[cell]
                        .iter()
                        .all(|&v| mesh.verts[v].into_iter().all(|x| x > 0.125 && x < 0.875))
                })
        })
        .unwrap();
    let velocity = counterflow(&system, edge, 0.02);
    let depth = 1.7;
    let plan = system
        .prepare_density_transport(&velocity, depth, 1e-12, 1 << 22, || false)
        .unwrap();
    let rho = |point: [f64; 2]| 1.0 + 8.0 * point[0] + 3.0 * point[1];
    let specific = |point: [f64; 2]| -20.0 + 3.0 * point[0] - point[1];
    let cells: Vec<_> = mesh
        .centroids
        .iter()
        .enumerate()
        .map(|(cell, &point)| {
            // The physical specific field is affine, but its cell state is
            // the integral of rho*s divided by mass, not s(geometric center).
            let mean = tri_quad(mesh.tris[cell].map(|v| mesh.verts[v]), mesh.areas[cell])
                .into_iter()
                .map(|(point, weight)| weight * rho(point) * specific(point))
                .sum::<f64>()
                / (mesh.areas[cell] * rho(point));
            Cell {
                density_kg_m3: rho(point),
                specific: [mean, 7.0],
            }
        })
        .collect();
    let mut work = plan.workspace::<2>(1 << 22).unwrap();
    let mut output = cells.clone();
    let dt = 0.01 * plan.maximum_timestep_s(0.6).unwrap();
    let report = plan
        .advance_into(
            &cells,
            &vec![None; mesh.edges.len()],
            dt,
            0.6,
            &mut output,
            &mut work,
            || false,
        )
        .unwrap();
    let face = &mesh.edges[edge];
    let a = mesh.verts[face.verts.0];
    let b = mesh.verts[face.verts.1];
    for (direction, donor, left, right) in [(0, face.tris.0, 0.5, 1.0), (1, face.tris.1, 0.0, 0.5)]
    {
        // Independent two-point Gauss exactly integrates the PHYSICAL cubic
        // rho*s*u_n on each counterflow interval. A geometric-center stencil
        // distorts the affine field and fails this oracle at variable density.
        let mut expected = 0.0;
        for t in [-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()] {
            let s = 0.5 * (left + right) + 0.5 * (right - left) * t;
            let point = core::array::from_fn(|k| a[k] + s * (b[k] - a[k]));
            let normal = (12.0 * velocity[2 * edge + 1] * (s - 0.5)).abs();
            expected += 0.5 * (right - left) * normal * rho(point) * specific(point);
        }
        let accepted = work.face_transfers().unwrap()[edge];
        near(
            accepted.carried[direction][0],
            dt * depth * face.len * expected,
            2e-14,
        );
        assert!(
            (accepted.carried[direction][0]
                - accepted.mass_kg[direction] * cells[donor].specific[0])
                .abs()
                > 1e-8
        );
    }
    near(report.relative_mass_balance_residual, 0.0, 3e-15);
    for residual in report.relative_carried_balance_residual {
        near(residual, 0.0, 3e-15);
    }
    assert!(output.iter().all(|state| state.specific[1] == 7.0));
}

#[test]
fn inventory_limiter_preserves_species_bounds_sums_and_mass_flux_at_density_contrast() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 10, 3));
    let system = FluxSystem::new(&mesh)
        .with_density_reconstruction(DensityReconstruction::BoundedLinearInventories);
    let mut velocity = uniform_velocity(&system, [0.5, 0.0]);
    // Add zero-mean counterflow on EVERY face, including physical boundaries.
    for edge in 0..mesh.edges.len() {
        velocity[2 * edge + 1] = 0.15;
    }
    let plan = system
        .prepare_density_transport(&velocity, 2.5, 1e-12, 1 << 22, || false)
        .unwrap();
    let cells: Vec<_> = mesh
        .centroids
        .iter()
        .map(|point| {
            let species = if point[0] < 0.5 { 1.0 } else { 0.0 };
            Cell {
                density_kg_m3: if point[0] < 0.5 { 1000.0 } else { 1e-6 },
                specific: [species, 1.0 - species, 1_099_511_627_776.0],
            }
        })
        .collect();
    let inflow: Vec<_> = mesh
        .edges
        .iter()
        .map(|edge| {
            (edge.tris.1 == usize::MAX).then_some(Cell {
                density_kg_m3: 1000.0,
                specific: [0.25, 0.75, 1_099_511_627_776.0],
            })
        })
        .collect();
    let mut output = cells.clone();
    let mut work = plan.ssprk2_workspace::<3>(1 << 22).unwrap();
    let dt = 0.1;
    let report = work
        .advance_into(
            &cells,
            &inflow,
            dt,
            DensitySspRk2Controls {
                maximum_courant: 0.9,
                maximum_substeps: 256,
            },
            &mut output,
            || false,
        )
        .unwrap()
        .transport;
    for state in &output {
        assert!(state.density_kg_m3 > 0.0);
        for value in &state.specific[..2] {
            assert!((-2e-13..=1.0 + 2e-13).contains(value));
        }
        near(state.specific[0] + state.specific[1], 1.0, 3e-13);
        assert_eq!(state.specific[2], 1_099_511_627_776.0);
    }
    near(report.relative_mass_balance_residual, 0.0, 3e-14);
    for residual in report.relative_carried_balance_residual {
        near(residual, 0.0, 5e-14);
    }
    let flux = work.prepare_mass_flux(1 << 22, || false).unwrap();
    for (edge, transfer) in work.face_transfers().unwrap().iter().enumerate() {
        for direction in 0..2 {
            near(
                dt * 2.5 * flux.edge_directional_mass_rates(edge).unwrap()[direction],
                transfer.mass_kg[direction],
                3e-13 * transfer.mass_kg[direction].max(1.0),
            );
            near(
                transfer.carried[direction][0] + transfer.carried[direction][1],
                transfer.mass_kg[direction],
                3e-13 * transfer.mass_kg[direction].max(1.0),
            );
        }
    }
}

#[test]
fn inventory_storage_and_late_cancellation_preserve_transactional_output() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 5, 2));
    let system = FluxSystem::new(&mesh)
        .with_density_reconstruction(DensityReconstruction::BoundedLinearInventories);
    let plan = system
        .prepare_density_transport(
            &uniform_velocity(&system, [0.5, 0.0]),
            1.0,
            1e-12,
            1 << 22,
            || false,
        )
        .unwrap();
    let cells: Vec<_> = (0..mesh.tris.len())
        .map(|cell| thermal_average(&mesh, cell, 0.0))
        .collect();
    let inflow: Vec<_> = mesh
        .edges
        .iter()
        .map(|edge| {
            (edge.tris.1 == usize::MAX).then_some(Cell {
                density_kg_m3: 1.0,
                specific: [200_000.0, 0.3, -123.0],
            })
        })
        .collect();
    let required = match plan.workspace::<3>(0) {
        Err(Error::Storage { required, .. }) => required,
        _ => panic!("bounded workspace must admit its complete storage"),
    };
    assert!(matches!(
        plan.workspace::<3>(required - 1),
        Err(Error::Storage { .. })
    ));
    let mut work = plan.workspace::<3>(required).unwrap();
    let mut output = cells.clone();
    let dt = 0.5 * plan.maximum_timestep_s(0.6).unwrap();
    let mut calls = 0;
    plan.advance_into(&cells, &inflow, dt, 0.6, &mut output, &mut work, || {
        calls += 1;
        false
    })
    .unwrap();
    let accepted = output.clone();
    let transfers = work.face_transfers().unwrap().to_vec();
    let mut polls = 0;
    assert_eq!(
        plan.advance_into(&cells, &inflow, dt, 0.6, &mut output, &mut work, || {
            polls += 1;
            polls == calls
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(output, accepted);
    assert!(work.face_transfers().is_none());
    plan.advance_into(&cells, &inflow, dt, 0.6, &mut output, &mut work, || false)
        .unwrap();
    assert_eq!(output, accepted);
    assert_eq!(work.face_transfers().unwrap(), transfers);
}
