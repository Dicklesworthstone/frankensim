//! G1/G3: independently derived SSPRK2 transport and accepted flux weights.

use super::{Cell, Error, FluxSystem, Mesh2, TriMesh, counterflow, near};
use fs_flux::ns::transport::density::ssprk2::DensitySspRk2Controls;

fn controls() -> DensitySspRk2Controls {
    DensitySspRk2Controls {
        maximum_courant: 0.3,
        maximum_substeps: 128,
    }
}

#[test]
fn closed_rk_stages_match_exact_recurrence_and_weighted_directional_flux() {
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
    let datum = 1_099_511_627_776.0;
    let mut old = [Cell {
        density_kg_m3: 2.0,
        specific: [3.0, datum],
    }; 2];
    old[b] = Cell {
        density_kg_m3: 8.0,
        specific: [-1.0, datum],
    };
    let inflow = vec![None; mesh.edges.len()];
    let mut next = old;
    let mut work = plan.ssprk2_workspace::<2>(1 << 20).unwrap();
    let result = work
        .advance_into(&old, &inflow, 5.0, controls(), &mut next, || false)
        .unwrap();
    assert_eq!(result.substeps, 8);
    assert_eq!(result.stage_evaluations, 16);
    assert_eq!(result.substep_timestep_s, 0.625);
    near(result.maximum_stage_courant, 0.25, 3e-15);
    near(result.transport.maximum_courant, 2.0, 3e-15);
    // Independent two-cell generator: one Euler difference multiplier is
    // 1-2*alpha=1/2. SSPRK2 instead has (1+(1/2)^2)/2=5/8.
    let r = (5.0_f64 / 8.0).powi(8);
    near(next[a].density_kg_m3, 5.0 - 3.0 * r, 5e-14);
    near(next[b].density_kg_m3, 5.0 + 3.0 * r, 5e-14);
    near(
        next[a].density_kg_m3 * next[a].specific[0],
        -1.0 + 7.0 * r,
        8e-14,
    );
    near(
        next[b].density_kg_m3 * next[b].specific[0],
        -1.0 - 7.0 * r,
        8e-14,
    );
    assert!(next.iter().all(|cell| cell.specific[1] == datum));
    let mean = work.average_donor_density().unwrap();
    near(mean[a].density_kg_m3, 17.0 / 4.0 + 3.0 * r / 4.0, 5e-14);
    near(mean[b].density_kg_m3, 23.0 / 4.0 - 3.0 * r / 4.0, 5e-14);
    let transfer = work.face_transfers().unwrap()[edge];
    near(transfer.mass_kg[0], 17.0 / 2.0 + 3.0 * r / 2.0, 1e-13);
    near(transfer.mass_kg[1], 23.0 / 2.0 - 3.0 * r / 2.0, 1e-13);
    near(transfer.carried[0][0], 3.0 / 2.0 - 7.0 * r / 2.0, 1e-13);
    near(transfer.carried[1][0], -11.0 / 2.0 + 7.0 * r / 2.0, 1e-13);
    near(result.transport.mass_after_kg, 10.0, 5e-14);
    near(result.transport.carried_after[0], -2.0, 1e-13);
    near(result.transport.relative_mass_balance_residual, 0.0, 6e-15);
    for residual in result.transport.relative_carried_balance_residual {
        near(residual, 0.0, 1e-14);
    }
    let flux = work.prepare_mass_flux(1 << 20, || false).unwrap();
    let rates = flux.edge_directional_mass_rates(edge).unwrap();
    for direction in 0..2 {
        near(10.0 * rates[direction], transfer.mass_kg[direction], 2e-13);
    }
    for cell in 0..2 {
        near(
            next[cell].density_kg_m3 - old[cell].density_kg_m3
                + 5.0 * flux.cell_divergence(cell).unwrap(),
            0.0,
            3e-13,
        );
    }
}

#[test]
fn fixed_spatial_density_and_extensive_inventory_converge_at_second_temporal_order() {
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
    let mut old = [Cell {
        density_kg_m3: 2.0,
        specific: [3.0],
    }; 2];
    old[b] = Cell {
        density_kg_m3: 8.0,
        specific: [-1.0],
    };
    let inflow = vec![None; mesh.edges.len()];
    let mut work = plan.ssprk2_workspace::<1>(1 << 20).unwrap();
    // Spatial discretization is FIXED: d(rho_a-rho_b)/dt=-0.8*(rho_a-rho_b),
    // and exactly the same ODE holds for the extensive E difference.
    let exact = (-0.8_f64).exp();
    let mut errors = Vec::new();
    for count in [4usize, 8, 16] {
        let mut next = old;
        let report = work
            .advance_into(
                &old,
                &inflow,
                1.0,
                DensitySspRk2Controls {
                    maximum_courant: 0.44 / count as f64,
                    maximum_substeps: count,
                },
                &mut next,
                || false,
            )
            .unwrap();
        assert_eq!(report.substeps, count);
        assert_eq!(report.stage_evaluations, 2 * count);
        errors.push([
            (next[a].density_kg_m3 - next[b].density_kg_m3 + 6.0 * exact).abs(),
            (next[a].density_kg_m3 * next[a].specific[0]
                - next[b].density_kg_m3 * next[b].specific[0]
                - 14.0 * exact)
                .abs(),
        ]);
    }
    for pair in errors.windows(2) {
        for component in 0..2 {
            let ratio = pair[0][component] / pair[1][component];
            assert!(
                (3.8..4.6).contains(&ratio),
                "ratio={ratio:e}, errors={errors:?}"
            );
        }
    }
    // SSP convexity also covers admitted stages beyond the Euler difference's
    // sign reversal at alpha=1/2. No artificial monotone-Euler restriction.
    let mut next = old;
    let result = work
        .advance_into(
            &old,
            &inflow,
            2.0,
            DensitySspRk2Controls {
                maximum_courant: 0.9,
                maximum_substeps: 1,
            },
            &mut next,
            || false,
        )
        .unwrap();
    assert_eq!(result.substeps, 1);
    near(next[a].density_kg_m3, 2.96, 4e-14);
    near(next[b].density_kg_m3, 7.04, 4e-14);
    assert!(
        next.iter()
            .all(|cell| cell.density_kg_m3 > 0.0 && (-1.0..=3.0).contains(&cell.specific[0]))
    );
}

#[test]
fn throughflow_uses_rk_boundary_weights_without_double_counting_stage_time() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let flow: Vec<_> = mesh
        .edges
        .iter()
        .flat_map(|edge| [edge.normal[0], 0.0])
        .collect();
    let plan = system
        .prepare_density_transport(&flow, 2.0, 1e-14, 1 << 20, || false)
        .unwrap();
    let a = usize::from(mesh.centroids[1][0] < mesh.centroids[0][0]);
    let b = 1 - a;
    let mut old = [Cell {
        density_kg_m3: 2.0,
        specific: [3.0],
    }; 2];
    old[b] = Cell {
        density_kg_m3: 8.0,
        specific: [-1.0],
    };
    let inflow: Vec<_> = mesh
        .edges
        .iter()
        .map(|edge| {
            (edge.tris.1 == usize::MAX).then_some(Cell {
                density_kg_m3: 10.0,
                specific: [5.0],
            })
        })
        .collect();
    let mut next = old;
    let mut work = plan.ssprk2_workspace::<1>(1 << 20).unwrap();
    let result = work
        .advance_into(&old, &inflow, 1.0, controls(), &mut next, || false)
        .unwrap();
    assert_eq!(result.substeps, 8);
    // Independent Jordan-block RK polynomial: relative to the fixed boundary
    // state its diagonal is r=25/32 and lower off-diagonal s=3/16.
    let r = 25.0_f64 / 32.0;
    let s = 3.0 / 16.0;
    let rho = [
        10.0 - 8.0 * r.powi(8),
        10.0 - 2.0 * r.powi(8) - 64.0 * s * r.powi(7),
    ];
    let extensive = [
        50.0 - 44.0 * r.powi(8),
        50.0 - 58.0 * r.powi(8) - 352.0 * s * r.powi(7),
    ];
    near(next[a].density_kg_m3, rho[0], 8e-14);
    near(next[b].density_kg_m3, rho[1], 8e-14);
    near(
        next[a].density_kg_m3 * next[a].specific[0],
        extensive[0],
        5e-13,
    );
    near(
        next[b].density_kg_m3 * next[b].specific[0],
        extensive[1],
        5e-13,
    );
    let report = result.transport;
    near(report.inlet_mass_kg, 20.0, 1e-13);
    near(report.carried_in[0], 100.0, 8e-13);
    near(report.outlet_mass_kg, 30.0 - rho[0] - rho[1], 2e-13);
    near(
        report.carried_out[0],
        98.0 - extensive[0] - extensive[1],
        1e-12,
    );
    near(report.relative_mass_balance_residual, 0.0, 1e-14);
    near(report.relative_carried_balance_residual[0], 0.0, 2e-14);
    let internal = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let direction = usize::from(mesh.edges[internal].tris.0 != a);
    let transfer = work.face_transfers().unwrap()[internal];
    near(transfer.mass_kg[direction], 22.0 - rho[0], 2e-13);
    near(transfer.carried[direction][0], 106.0 - extensive[0], 1e-12);
    let flux = work.prepare_mass_flux(1 << 20, || false).unwrap();
    for cell in 0..2 {
        near(
            next[cell].density_kg_m3 - old[cell].density_kg_m3
                + flux.cell_divergence(cell).unwrap(),
            0.0,
            3e-13,
        );
    }
    for edge in 0..mesh.edges.len() {
        let rates = flux.edge_directional_mass_rates(edge).unwrap();
        for direction in 0..2 {
            near(
                2.0 * rates[direction],
                work.face_transfers().unwrap()[edge].mass_kg[direction],
                2e-13,
            );
        }
    }
}

#[test]
fn admitted_extreme_density_stages_preserve_a_light_donors_finite_inventory() {
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
    let mut old = [Cell {
        density_kg_m3: 1e-300,
        specific: [1e100, -1e100],
    }; 2];
    old[b] = Cell {
        density_kg_m3: 1e100,
        specific: [0.0; 2],
    };
    assert_eq!(old[a].density_kg_m3 / old[b].density_kg_m3, 0.0);
    let mut next = old;
    let mut work = plan.ssprk2_workspace::<2>(1 << 20).unwrap();
    let result = work
        .advance_into(
            &old,
            &vec![None; mesh.edges.len()],
            0.625,
            controls(),
            &mut next,
            || false,
        )
        .unwrap();
    assert_eq!(result.substeps, 1);
    // Accepted RK mixing fractions are 13/16 and 3/16 at alpha=1/4.
    for (cell, fraction) in [(a, 13.0 / 16.0), (b, 3.0 / 16.0)] {
        for (component, sign) in [(0, 1.0), (1, -1.0)] {
            near(
                next[cell].density_kg_m3 * next[cell].specific[component] / (sign * 1e-200),
                fraction,
                2e-15,
            );
        }
    }
    for (component, sign) in [(0, 1.0), (1, -1.0)] {
        near(
            result.transport.carried_after[component] / (sign * 1e-200),
            1.0,
            3e-15,
        );
        near(
            result.transport.relative_carried_balance_residual[component],
            0.0,
            3e-15,
        );
    }
}

#[test]
fn equal_subnormal_stage_face_amounts_survive_the_rk_half_weights() {
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
    let tiny = f64::from_bits(1);
    let old = [Cell {
        density_kg_m3: 4.0 * tiny,
        specific: [1.0],
    }; 2];
    let mut next = old;
    let mut work = plan.ssprk2_workspace::<1>(1 << 20).unwrap();
    let result = work
        .advance_into(
            &old,
            &vec![None; mesh.edges.len()],
            0.625,
            controls(),
            &mut next,
            || false,
        )
        .unwrap();
    assert_eq!(result.substeps, 1);
    assert_eq!(next, old);
    let face = work.face_transfers().unwrap()[edge];
    assert_eq!(face.mass_kg, [tiny; 2]);
    assert_eq!(face.carried, [[tiny]; 2]);
    for mean in work.average_donor_density().unwrap() {
        assert_eq!(mean.density_kg_m3, 4.0 * tiny);
    }
    assert_eq!(result.transport.mass_balance_residual_kg, 0.0);
}

#[test]
fn bounded_rk_work_and_late_cancellation_are_atomic_and_retryable() {
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
        plan.ssprk2_workspace::<3>(0),
        Err(Error::Storage { .. })
    ));
    let old = super::initial();
    let inflow = vec![None; mesh.edges.len()];
    let mut output = old;
    let mut work = plan.ssprk2_workspace::<3>(1 << 20).unwrap();
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
    assert_eq!(retried, accepted);
    assert_eq!(output, expected);
    assert_eq!(work.face_transfers().unwrap(), faces);
    output = sentinel;
    assert!(matches!(
        work.advance_into(
            &old,
            &inflow,
            5.0,
            DensitySspRk2Controls {
                maximum_substeps: 4,
                ..controls()
            },
            &mut output,
            || false,
        ),
        Err(Error::SubstepLimit { maximum: 4 })
    ));
    assert_eq!(output, sentinel);
    assert!(work.face_transfers().is_none());
}
