#![cfg(feature = "continuum")]
//! G1/G3: original BDM face streams drive conservative changing cell density,
//! momentum and enthalpy; no independent net-flow or fixed-mass surrogate.

use fs_flux::ns::transport::density::{
    DensityTransportCell as Cell, DensityTransportError as Error,
};
use fs_flux::{FluxSystem, TriMesh};
use fs_solid::Mesh2;

#[path = "density_transport/mass_flux.rs"]
mod mass_flux;

#[path = "density_transport/momentum.rs"]
mod momentum;

#[path = "density_transport/subcycle.rs"]
mod subcycle;

#[path = "density_transport/remap.rs"]
mod remap;

#[path = "density_transport/ssprk2.rs"]
mod ssprk2;

#[path = "density_transport/linear.rs"]
mod linear;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

fn counterflow(system: &FluxSystem<'_>, edge: usize, rate_m2_s: f64) -> Vec<f64> {
    let mut velocity = vec![0.0; system.n_u];
    // m0=0, u.n=12*m1*(s-1/2): each lobe integrates to 1.5*m1*length.
    velocity[2 * edge + 1] = rate_m2_s / (1.5 * system.mesh.edges[edge].len);
    velocity
}

fn initial() -> [Cell<3>; 2] {
    [
        Cell {
            density_kg_m3: 2.0,
            specific: [3.0, -2.0, 10.0],
        },
        Cell {
            density_kg_m3: 8.0,
            specific: [-1.0, 4.0, 100.0],
        },
    ]
}

#[test]
fn counterflow_transports_unequal_mass_momentum_and_enthalpy_from_one_face() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let edge = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let velocity = counterflow(&system, edge, 0.2);
    let plan = system
        .prepare_density_transport(&velocity, 2.0, 1e-14, 1 << 20, || false)
        .unwrap();
    assert_eq!(plan.cell_volumes_m3(), &[1.0, 1.0]);
    assert_eq!(
        plan.edge_volume_rates(),
        system.edge_volume_transport(&velocity).unwrap()
    );
    assert_eq!(plan.edge_volume_rates()[edge].net_m2_s(), 0.0);
    let old = initial();
    let mut output = [Cell::default(); 2];
    let mut scratch = plan.workspace::<3>(1 << 20).unwrap();
    assert!(scratch.face_transfers().is_none());
    let report = plan
        .advance_into(
            &old,
            &vec![None; mesh.edges.len()],
            0.5,
            0.8,
            &mut output,
            &mut scratch,
            || false,
        )
        .unwrap();
    let (i, j) = mesh.edges[edge].tris;
    let transfers = scratch.face_transfers().unwrap();
    for (direction, donor) in [(0, i), (1, j)] {
        let mass = 0.2 * old[donor].density_kg_m3;
        near(transfers[edge].mass_kg[direction], mass, 4e-15);
        for k in 0..3 {
            near(
                transfers[edge].carried[direction][k],
                mass * old[donor].specific[k],
                3e-13,
            );
        }
    }
    for (cell, other) in [(i, j), (j, i)] {
        let mass = 0.8 * old[cell].density_kg_m3 + 0.2 * old[other].density_kg_m3;
        near(output[cell].density_kg_m3, mass, 4e-15);
        for k in 0..3 {
            let inventory = 0.8 * old[cell].density_kg_m3 * old[cell].specific[k]
                + 0.2 * old[other].density_kg_m3 * old[other].specific[k];
            near(mass * output[cell].specific[k], inventory, 2e-13);
        }
    }
    near(report.mass_before_kg, 10.0, 0.0);
    near(report.mass_after_kg, 10.0, 2e-15);
    for (k, total) in [-2.0, 28.0, 820.0].into_iter().enumerate() {
        near(report.carried_before[k], total, 0.0);
        near(report.carried_after[k], total, 2e-13);
        near(report.carried_balance_residual[k], 0.0, 2e-13);
        near(report.relative_carried_balance_residual[k], 0.0, 5e-16);
    }
    let kinetic = |states: &[Cell<3>]| {
        states
            .iter()
            .map(|x| 0.5 * x.density_kg_m3 * (x.specific[0].powi(2) + x.specific[1].powi(2)))
            .sum::<f64>()
    };
    assert!(kinetic(&output) < kinetic(&old));
    // Numerical momentum mixing is not silently added to thermal enthalpy.
    near(report.carried_after[2], report.carried_before[2], 2e-13);
    assert_eq!(old, initial());
}

#[test]
fn high_density_contrast_does_not_cancel_the_small_donors_finite_inventory() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let edge = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let plan = system
        .prepare_density_transport(
            &counterflow(&system, edge, 0.2),
            2.0,
            1e-14,
            1 << 20,
            || false,
        )
        .unwrap();
    let initial = [
        Cell {
            density_kg_m3: 1.0,
            specific: [1e30],
        },
        Cell {
            density_kg_m3: 1e30,
            specific: [0.0],
        },
    ];
    let mut output = initial;
    let mut scratch = plan.workspace::<1>(1 << 20).unwrap();
    let report = plan
        .advance_into(
            &initial,
            &vec![None; mesh.edges.len()],
            0.5,
            0.8,
            &mut output,
            &mut scratch,
            || false,
        )
        .unwrap();
    // In cell 0, 0.8*1*1e30 is retained beside 0.2*1e30 incoming mass
    // with zero carried quantity. A naive old + delta/m update returns zero.
    near(output[0].specific[0], 4.0, 3e-15);
    near(output[1].specific[0], 0.25, 3e-16);
    near(report.carried_after[0] / 1e30, 1.0, 4e-16);
    near(report.relative_carried_balance_residual[0], 0.0, 8e-16);
}

#[test]
fn density_and_extensive_enthalpy_converge_to_independent_exchange_solution() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let edge = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    // Q=depth*0.25=0.5 m³/s, V=1 m³. The difference decays as exp(-t).
    let plan = system
        .prepare_density_transport(
            &counterflow(&system, edge, 0.25),
            2.0,
            1e-14,
            1 << 20,
            || false,
        )
        .unwrap();
    let inflow = vec![None; mesh.edges.len()];
    let mut errors = Vec::new();
    for steps in [8, 16, 32, 64] {
        let mut current = initial();
        let mut next = current;
        let mut scratch = plan.workspace::<3>(1 << 20).unwrap();
        for _ in 0..steps {
            let report = plan
                .advance_into(
                    &current,
                    &inflow,
                    1.0 / f64::from(steps),
                    0.8,
                    &mut next,
                    &mut scratch,
                    || false,
                )
                .unwrap();
            near(report.relative_mass_balance_residual, 0.0, 4e-16);
            for residual in report.relative_carried_balance_residual {
                near(residual, 0.0, 8e-16);
            }
            core::mem::swap(&mut current, &mut next);
        }
        let exact_rho = 5.0 - 3.0 * (-1.0_f64).exp();
        let exact_enthalpy = 410.0 - 390.0 * (-1.0_f64).exp();
        errors.push([
            (current[0].density_kg_m3 - exact_rho).abs(),
            (current[0].density_kg_m3 * current[0].specific[2] - exact_enthalpy).abs(),
        ]);
        near(current.iter().map(|x| x.density_kg_m3).sum(), 10.0, 3e-14);
        for cell in current {
            assert!((2.0..=8.0).contains(&cell.density_kg_m3));
        }
    }
    for pair in errors.windows(2) {
        for k in 0..2 {
            assert!(pair[0][k] / pair[1][k] > 1.95);
        }
    }
}

#[test]
fn native_throughflow_preserves_specific_free_stream_as_density_changes() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(2.0, 1.0, 4, 2));
    let system = FluxSystem::new(&mesh);
    let velocity: Vec<_> = mesh.edges.iter().flat_map(|e| [e.normal[0], 0.0]).collect();
    let plan = system
        .prepare_density_transport(&velocity, 0.3, 1e-14, 1 << 20, || false)
        .unwrap();
    let specific = [7.0, -3.0, 1_099_511_627_776.0];
    let mut current: Vec<_> = mesh
        .centroids
        .iter()
        .map(|p| Cell {
            density_kg_m3: if p[0] < 1.0 { 0.001 } else { 1000.0 },
            specific,
        })
        .collect();
    let original = current.clone();
    let inflow: Vec<_> = mesh
        .edges
        .iter()
        .map(|e| {
            (e.tris.1 == usize::MAX).then_some(Cell {
                density_kg_m3: 3.0,
                specific,
            })
        })
        .collect();
    let mut next = current.clone();
    let mut scratch = plan.workspace::<3>(1 << 20).unwrap();
    let dt = 0.99 * plan.maximum_timestep_s(0.8).unwrap();
    for _ in 0..30 {
        let report = plan
            .advance_into(&current, &inflow, dt, 0.8, &mut next, &mut scratch, || {
                false
            })
            .unwrap();
        near(report.relative_mass_balance_residual, 0.0, 2e-15);
        for c in &next {
            assert_eq!(c.specific, specific);
            assert!(c.density_kg_m3 >= 0.001 * (1.0 - 1e-13));
            assert!(c.density_kg_m3 <= 1000.0 * (1.0 + 1e-13));
        }
        core::mem::swap(&mut current, &mut next);
    }
    assert_ne!(current, original);
}

#[test]
fn reversing_open_face_uses_both_actual_donors_and_closes_boundary_budgets() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let edge = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 == usize::MAX && e.tris.0 == 0)
        .unwrap();
    let plan = system
        .prepare_density_transport(
            &counterflow(&system, edge, 0.1),
            2.0,
            1e-14,
            1 << 20,
            || false,
        )
        .unwrap();
    let mut inflow = vec![None; mesh.edges.len()];
    inflow[edge] = Some(Cell {
        density_kg_m3: 10.0,
        specific: [10.0, 0.0, 200.0],
    });
    let old = initial();
    let mut next = old;
    let mut scratch = plan.workspace::<3>(1 << 20).unwrap();
    let report = plan
        .advance_into(&old, &inflow, 0.5, 0.8, &mut next, &mut scratch, || false)
        .unwrap();
    near(report.inlet_mass_kg, 1.0, 1e-15);
    near(report.outlet_mass_kg, 0.2, 2e-16);
    near(next[0].density_kg_m3, 2.8, 1e-15);
    assert_eq!(next[1], old[1]);
    for k in 0..3 {
        let expected =
            0.9 * old[0].density_kg_m3 * old[0].specific[k] + inflow[edge].unwrap().specific[k];
        near(next[0].density_kg_m3 * next[0].specific[k], expected, 2e-13);
        near(
            report.carried_after[k] - report.carried_before[k],
            report.carried_in[k] - report.carried_out[k],
            3e-13,
        );
    }
    let flux = scratch.face_transfers().unwrap()[edge];
    near(flux.mass_kg[0], 0.2, 2e-16);
    near(flux.mass_kg[1], 1.0, 1e-15);
    near(report.mass_balance_residual_kg, 0.0, 4e-16);
}

#[test]
fn physical_depth_density_and_specific_datum_scaling_are_coherent() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let edge = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let velocity = counterflow(&system, edge, 0.2);
    let mut results = Vec::new();
    for (depth, density_scale, datum) in [
        (2.0, 1.0, 0.0),
        (14.0, 1.0, 0.0),
        (2.0, 5.0, 0.0),
        (2.0, 1.0, 1_099_511_627_776.0),
    ] {
        let plan = system
            .prepare_density_transport(&velocity, depth, 1e-14, 1 << 20, || false)
            .unwrap();
        let states = initial().map(|mut cell| {
            cell.density_kg_m3 *= density_scale;
            cell.specific[2] += datum;
            cell
        });
        let mut next = states;
        let mut scratch = plan.workspace::<3>(1 << 20).unwrap();
        let report = plan
            .advance_into(
                &states,
                &vec![None; mesh.edges.len()],
                0.5,
                0.8,
                &mut next,
                &mut scratch,
                || false,
            )
            .unwrap();
        for cell in &mut next {
            cell.density_kg_m3 /= density_scale;
            cell.specific[2] -= datum;
        }
        results.push((next, report));
    }
    for (next, _) in &results[1..] {
        for c in 0..2 {
            near(next[c].density_kg_m3, results[0].0[c].density_kg_m3, 2e-15);
            for k in 0..3 {
                near(
                    next[c].specific[k],
                    results[0].0[c].specific[k],
                    if k == 2 { 2e-4 } else { 2e-15 },
                );
            }
        }
    }
    near(
        results[1].1.mass_after_kg,
        7.0 * results[0].1.mass_after_kg,
        3e-14,
    );
    near(
        results[2].1.mass_after_kg,
        5.0 * results[0].1.mass_after_kg,
        2e-14,
    );
}

#[test]
fn refusals_preserve_output_and_invalidate_any_previous_accepted_transfer_view() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let edge = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let mut velocity = counterflow(&system, edge, 0.2);
    assert!(matches!(
        system.prepare_density_transport(&velocity, 2.0, 1e-14, 0, || false),
        Err(Error::Storage { .. })
    ));
    let plan = system
        .prepare_density_transport(&velocity, 2.0, 1e-14, 1 << 20, || false)
        .unwrap();
    assert!(matches!(plan.workspace::<3>(0), Err(Error::Storage { .. })));
    let mut scratch = plan.workspace::<3>(1 << 20).unwrap();
    let old = initial();
    let inflow = vec![None; mesh.edges.len()];
    let mut output = old;
    plan.advance_into(&old, &inflow, 0.1, 0.8, &mut output, &mut scratch, || false)
        .unwrap();
    assert!(scratch.face_transfers().is_some());
    let sentinel = output;
    assert!(matches!(
        plan.advance_into(&old, &inflow, 100.0, 0.8, &mut output, &mut scratch, || {
            false
        }),
        Err(Error::Courant { .. })
    ));
    assert_eq!(output, sentinel);
    assert!(scratch.face_transfers().is_none());
    let mut count = 0;
    plan.advance_into(&old, &inflow, 0.1, 0.8, &mut output, &mut scratch, || {
        count += 1;
        false
    })
    .unwrap();
    let sentinel = output;
    let mut calls = 0;
    assert_eq!(
        plan.advance_into(&old, &inflow, 0.1, 0.8, &mut output, &mut scratch, || {
            calls += 1;
            calls == count
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(output, sentinel);
    assert!(scratch.face_transfers().is_none());
    let mut bad = old;
    bad[1].density_kg_m3 = -1.0;
    assert!(
        plan.advance_into(&bad, &inflow, 0.1, 0.8, &mut output, &mut scratch, || false)
            .is_err()
    );
    assert_eq!(output, sentinel);
    // A nonzero one-way interior flow does not conserve either cell's volume.
    velocity[2 * edge] = 1.0;
    assert!(matches!(
        system.prepare_density_transport(&velocity, 2.0, 1e-14, 1 << 20, || false),
        Err(Error::VolumeDivergence { .. })
    ));
}

#[test]
fn an_extra_global_counterflow_edge_cannot_bypass_the_original_cell_faces() {
    let mut mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let interior = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 != usize::MAX)
        .unwrap();
    let duplicate = mesh.edges.len();
    mesh.edges.push(mesh.edges[interior].clone());
    let system = FluxSystem::new(&mesh);
    // Each adjacent cell still sees zero net volume, but neither original
    // three-face stencil lists this extra global edge. It must not transport
    // mass which the extensive-mixture fallback would then omit.
    let velocity = counterflow(&system, duplicate, 0.2);
    assert!(matches!(
        system.prepare_density_transport(&velocity, 2.0, 1e-14, 1 << 20, || false),
        Err(Error::Invalid(
            "every global edge must belong to its original oriented cell faces"
        ))
    ));
}

#[test]
fn mass_only_rest_and_missing_reverse_boundary_state_are_explicit() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let still = system
        .prepare_density_transport(&vec![0.0; system.n_u], 2.0, 1e-14, 1 << 20, || false)
        .unwrap();
    assert_eq!(still.maximum_timestep_s(0.8).unwrap(), f64::INFINITY);
    let old = [
        Cell {
            density_kg_m3: 2.0,
            specific: [],
        },
        Cell {
            density_kg_m3: 8.0,
            specific: [],
        },
    ];
    let mut next = old;
    let mut scratch = still.workspace::<0>(1 << 20).unwrap();
    let boundary = vec![None; mesh.edges.len()];
    let report = still
        .advance_into(&old, &boundary, 100.0, 0.8, &mut next, &mut scratch, || {
            false
        })
        .unwrap();
    assert_eq!(next, old);
    assert_eq!(report.mass_balance_residual_kg, 0.0);
    let edge = mesh
        .edges
        .iter()
        .position(|e| e.tris.1 == usize::MAX)
        .unwrap();
    let open = system
        .prepare_density_transport(
            &counterflow(&system, edge, 0.1),
            2.0,
            1e-14,
            1 << 20,
            || false,
        )
        .unwrap();
    assert!(
        open.advance_into(&old, &boundary, 0.1, 0.8, &mut next, &mut scratch, || false)
            .is_err()
    );
    assert_eq!(next, old);
    assert!(scratch.face_transfers().is_none());
}
