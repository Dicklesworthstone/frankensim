//! G1/G3: exact donor-cell recurrences and one compatible macro mass flux.

use super::{Cell, Error, FluxSystem, Mesh2, TriMesh, counterflow, near};
use fs_flux::ns::transport::density::subcycle::DensitySubcycleControls;

fn controls() -> DensitySubcycleControls {
    DensitySubcycleControls {
        maximum_courant: 0.3,
        maximum_substeps: 128,
    }
}

#[test]
fn closed_subcycles_match_independent_exchange_and_integrated_flux() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let edge = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let (a, b) = mesh.edges[edge].tris;
    let flow = counterflow(&system, edge, 0.2);
    let plan = system
        .prepare_density_transport(&flow, 2.0, 1e-14, 1 << 20, || false)
        .unwrap();
    let mut old = [Cell::<1>::default(); 2];
    old[a] = Cell {
        density_kg_m3: 2.0,
        specific: [3.0],
    };
    old[b] = Cell {
        density_kg_m3: 8.0,
        specific: [-1.0],
    };
    let inflow = vec![None; mesh.edges.len()];
    let mut next = old;
    let mut work = plan.subcycle_workspace::<1>(1 << 20).unwrap();
    let result = work
        .advance_into(&old, &inflow, 5.0, controls(), &mut next, || false)
        .unwrap();
    assert_eq!(result.substeps, 8);
    assert_eq!(result.substep_timestep_s, 0.625);
    near(result.maximum_substep_courant, 0.25, 2e-15);
    near(result.transport.maximum_courant, 2.0, 2e-15);
    // Independently diagonalize [[3/4,1/4],[1/4,3/4]]. Its difference
    // eigenvalue is 1/2; both density and E=rho*s obey this same map.
    near(next[a].density_kg_m3, 1277.0 / 256.0, 4e-14);
    near(next[b].density_kg_m3, 1283.0 / 256.0, 4e-14);
    near(
        next[a].density_kg_m3 * next[a].specific[0],
        -249.0 / 256.0,
        4e-14,
    );
    near(
        next[b].density_kg_m3 * next[b].specific[0],
        -263.0 / 256.0,
        4e-14,
    );
    let mean = work.average_donor_density().unwrap();
    near(mean[a].density_kg_m3, 4355.0 / 1024.0, 3e-14);
    near(mean[b].density_kg_m3, 5885.0 / 1024.0, 3e-14);
    let actual = work.face_transfers().unwrap()[edge];
    near(actual.mass_kg[0], 4355.0 / 512.0, 8e-14);
    near(actual.mass_kg[1], 5885.0 / 512.0, 8e-14);
    near(actual.carried[0][0], 761.0 / 512.0, 8e-14);
    near(actual.carried[1][0], -2809.0 / 512.0, 8e-14);
    near(result.transport.mass_after_kg, 10.0, 6e-14);
    near(result.transport.carried_after[0], -2.0, 8e-14);
    near(result.transport.mass_balance_residual_kg, 0.0, 6e-14);
    near(result.transport.carried_balance_residual[0], 0.0, 8e-14);
    let flux = work.prepare_mass_flux(1 << 20, || false).unwrap();
    let directional = flux.edge_directional_mass_rates(edge).unwrap();
    for d in 0..2 {
        near(directional[d] * 5.0 * 2.0, actual.mass_kg[d], 1e-13);
    }
    for c in 0..2 {
        near(
            next[c].density_kg_m3 - old[c].density_kg_m3 + 5.0 * flux.cell_divergence(c).unwrap(),
            0.0,
            2e-13,
        );
    }
}

#[test]
fn throughflow_subcycles_keep_actual_boundary_mass_and_carried_exchange() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let flow: Vec<_> = mesh.edges.iter().flat_map(|e| [e.normal[0], 0.0]).collect();
    let plan = system
        .prepare_density_transport(&flow, 2.0, 1e-14, 1 << 20, || false)
        .unwrap();
    let upstream = usize::from(mesh.centroids[1][0] < mesh.centroids[0][0]);
    let downstream = 1 - upstream;
    let datum = 1_099_511_627_776.0;
    let mut old = [Cell::<2>::default(); 2];
    old[upstream] = Cell {
        density_kg_m3: 2.0,
        specific: [3.0, datum],
    };
    old[downstream] = Cell {
        density_kg_m3: 8.0,
        specific: [-1.0, datum],
    };
    let inflow: Vec<_> = mesh
        .edges
        .iter()
        .map(|e| {
            (e.tris.1 == usize::MAX).then_some(Cell {
                density_kg_m3: 10.0,
                specific: [5.0, datum],
            })
        })
        .collect();
    let mut next = old;
    let mut work = plan.subcycle_workspace::<2>(1 << 20).unwrap();
    let result = work
        .advance_into(&old, &inflow, 1.0, controls(), &mut next, || false)
        .unwrap();
    assert_eq!(result.substeps, 8);
    // Independent Jordan-block recurrence: x1'=3x1/4+xB/4,
    // x2'=3x2/4+x1/4, with xB=10 for rho and xB=50 for rho*s.
    near(next[upstream].density_kg_m3, 75359.0 / 8192.0, 6e-14);
    near(next[downstream].density_kg_m3, 251135.0 / 32768.0, 6e-14);
    near(
        next[upstream].density_kg_m3 * next[upstream].specific[0],
        747029.0 / 16384.0,
        4e-13,
    );
    near(
        next[downstream].density_kg_m3 * next[downstream].specific[0],
        1063219.0 / 32768.0,
        4e-13,
    );
    for cell in next {
        assert_eq!(cell.specific[1], datum);
    }
    let mean = work.average_donor_density().unwrap();
    near(mean[upstream].density_kg_m3, 104865.0 / 16384.0, 6e-14);
    near(mean[downstream].density_kg_m3, 430469.0 / 65536.0, 6e-14);
    let report = result.transport;
    near(report.inlet_mass_kg, 20.0, 1e-13);
    near(report.outlet_mass_kg, 430469.0 / 32768.0, 1e-13);
    near(report.carried_in[0], 100.0, 1e-12);
    near(report.carried_out[0], 653987.0 / 32768.0, 3e-13);
    near(report.relative_mass_balance_residual, 0.0, 3e-15);
    for residual in report.relative_carried_balance_residual {
        near(residual, 0.0, 5e-15);
    }
    let internal = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let direction = usize::from(mesh.edges[internal].tris.0 != upstream);
    let transfer = work.face_transfers().unwrap()[internal];
    near(transfer.mass_kg[direction], 104865.0 / 8192.0, 1e-13);
    near(transfer.carried[direction][0], 989675.0 / 16384.0, 6e-13);
    let flux = work.prepare_mass_flux(1 << 20, || false).unwrap();
    for cell in 0..2 {
        near(
            next[cell].density_kg_m3 - old[cell].density_kg_m3
                + flux.cell_divergence(cell).unwrap(),
            0.0,
            2e-13,
        );
    }
    for edge in 0..mesh.edges.len() {
        let mass = flux.edge_directional_mass_rates(edge).unwrap();
        for direction in 0..2 {
            near(
                2.0 * mass[direction],
                work.face_transfers().unwrap()[edge].mass_kg[direction],
                1e-13,
            );
        }
    }
}

#[test]
fn a_single_substep_preserves_original_euler_states_faces_and_flux() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let edge = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let flow = counterflow(&system, edge, 0.2);
    let plan = system
        .prepare_density_transport(&flow, 0.7, 1e-14, 1 << 20, || false)
        .unwrap();
    let old = super::initial();
    let inflow = vec![None; mesh.edges.len()];
    let mut direct = old;
    let mut next = old;
    let mut original = plan.workspace::<3>(1 << 20).unwrap();
    let expected = plan
        .advance_into(&old, &inflow, 0.1, 0.3, &mut direct, &mut original, || {
            false
        })
        .unwrap();
    let mut work = plan.subcycle_workspace::<3>(1 << 20).unwrap();
    let result = work
        .advance_into(&old, &inflow, 0.1, controls(), &mut next, || false)
        .unwrap();
    assert_eq!(result.substeps, 1);
    assert_eq!(next, direct);
    assert_eq!(work.face_transfers(), original.face_transfers());
    assert_eq!(result.transport.mass_before_kg, expected.mass_before_kg);
    assert_eq!(result.transport.mass_after_kg, expected.mass_after_kg);
    let a = work.prepare_mass_flux(1 << 20, || false).unwrap();
    let b = plan
        .prepare_mass_flux(&old, &inflow, 1 << 20, || false)
        .unwrap();
    for (cell, &point) in mesh.centroids.iter().enumerate() {
        assert_eq!(
            a.cell_value(cell, point).unwrap(),
            b.cell_value(cell, point).unwrap()
        );
    }
}

#[test]
fn scaled_density_mean_survives_large_and_subnormal_positive_inventories() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let edge = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let flow = counterflow(&system, edge, 0.5);
    // The physical masses stay representable. Raw sum(rho) overflows in the
    // heavy case, while sum(rho/n) would erase the light case completely.
    for (rho, depth) in [(1e308, 1e-280), (f64::from_bits(1), 1e300)] {
        let plan = system
            .prepare_density_transport(&flow, depth, 1e-14, 1 << 20, || false)
            .unwrap();
        let old = [Cell {
            density_kg_m3: rho,
            specific: [7.0],
        }; 2];
        let mut next = old;
        let mut work = plan.subcycle_workspace::<1>(1 << 20).unwrap();
        let result = work
            .advance_into(
                &old,
                &vec![None; mesh.edges.len()],
                4.0,
                DensitySubcycleControls {
                    maximum_courant: 0.75,
                    maximum_substeps: 8,
                },
                &mut next,
                || false,
            )
            .unwrap();
        assert_eq!(result.substeps, 8);
        for mean in work.average_donor_density().unwrap() {
            assert!(mean.density_kg_m3.is_finite() && mean.density_kg_m3 > 0.0);
            near(mean.density_kg_m3 / rho, 1.0, 3e-15);
        }
        for cell in next {
            assert_eq!(cell.specific, [7.0]);
        }
        near(result.transport.relative_mass_balance_residual, 0.0, 5e-15);
    }
}

#[test]
fn subcycle_budget_and_late_cancellation_preserve_output_and_allow_retry() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let edge = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let flow = counterflow(&system, edge, 0.2);
    let plan = system
        .prepare_density_transport(&flow, 2.0, 1e-14, 1 << 20, || false)
        .unwrap();
    assert!(matches!(
        plan.subcycle_workspace::<3>(0),
        Err(Error::Storage { .. })
    ));
    let old = super::initial();
    let inflow = vec![None; mesh.edges.len()];
    let mut output = old;
    let mut work = plan.subcycle_workspace::<3>(1 << 20).unwrap();
    let mut polls = 0;
    let accepted = work
        .advance_into(&old, &inflow, 5.0, controls(), &mut output, || {
            polls += 1;
            false
        })
        .unwrap();
    let expected = output;
    let faces = work.face_transfers().unwrap().to_vec();
    let sentinel = [Cell {
        density_kg_m3: 123.0,
        specific: [-4.0; 3],
    }; 2];
    output = sentinel;
    let mut seen = 0;
    assert!(matches!(
        work.advance_into(&old, &inflow, 5.0, controls(), &mut output, || {
            seen += 1;
            seen == polls
        }),
        Err(Error::Cancelled)
    ));
    assert_eq!(output, sentinel);
    assert!(work.face_transfers().is_none() && work.average_donor_density().is_none());
    assert!(matches!(
        work.prepare_mass_flux(1 << 20, || false),
        Err(Error::Invalid(_))
    ));
    let retried = work
        .advance_into(&old, &inflow, 5.0, controls(), &mut output, || false)
        .unwrap();
    assert_eq!(output, expected);
    assert_eq!(retried, accepted);
    assert_eq!(work.face_transfers().unwrap(), faces);
    output = sentinel;
    assert!(matches!(
        work.advance_into(
            &old,
            &inflow,
            5.0,
            DensitySubcycleControls {
                maximum_substeps: 4,
                ..controls()
            },
            &mut output,
            || false
        ),
        Err(Error::SubstepLimit { maximum: 4 })
    ));
    assert_eq!(output, sentinel);
    assert!(work.face_transfers().is_none());
    let mut invalid = old;
    invalid[0].density_kg_m3 = f64::NAN;
    assert!(matches!(
        work.advance_into(&invalid, &inflow, 5.0, controls(), &mut output, || false),
        Err(Error::Invalid(_))
    ));
    assert_eq!(output, sentinel);
}
