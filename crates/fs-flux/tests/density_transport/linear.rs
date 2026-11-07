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
