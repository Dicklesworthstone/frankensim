#![cfg(feature = "continuum")]
//! Physical G1/G3: density, native momentum and pressure from one IMEX step.

use fs_flux::bdm::{eval_basis, tri_quad};
use fs_flux::ns::material::mixed::density_transport::{
    ConservativeMomentumControls as Controls, ConservativeMomentumError as Error,
    ConservativeMomentumHistory as History,
};
use fs_flux::ns::material::{MaterialFluxNumerics, MaterialStokesBoundary as Boundary};
use fs_flux::ns::transport::density::{DensityTransportCell as Cell, DensityTransportError};
use fs_flux::{FluxSystem, TriMesh};
use fs_solid::Mesh2;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

fn controls(reference_density_kg_m3: f64) -> Controls {
    Controls {
        timestep_s: 0.01,
        depth_m: 0.7,
        reference_density_kg_m3,
        maximum_courant: 0.8,
        maximum_transport_storage_bytes: 1 << 20,
    }
}

fn constant(mesh: &TriMesh, velocity: [f64; 2]) -> Vec<f64> {
    let mut field = vec![0.0; 2 * mesh.edges.len()];
    for (index, edge) in mesh.edges.iter().enumerate() {
        field[2 * index] = velocity[0] * edge.normal[0] + velocity[1] * edge.normal[1];
    }
    field
}

fn momentum(system: &FluxSystem<'_>, density: &[Cell<0>], velocity: &[f64]) -> [f64; 2] {
    let mut momentum = [0.0; 2];
    for cell in 0..system.mesh.tris.len() {
        let p = system.mesh.tris[cell].map(|node| system.mesh.verts[node]);
        for (point, weight) in tri_quad(p, system.mesh.areas[cell]) {
            for local in 0..6 {
                let basis = eval_basis(&system.bases[cell], local, point);
                let coefficient =
                    velocity[2 * system.mesh.tri_edges[cell][local / 2].0 + local % 2];
                for axis in 0..2 {
                    momentum[axis] +=
                        density[cell].density_kg_m3 * weight * coefficient * basis[axis];
                }
            }
        }
    }
    momentum
}

#[test]
fn changing_density_preserves_uniform_flow_and_physical_reference_invariance() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 0.6, 3, 2));
    let system = FluxSystem::new(&mesh);
    let density: Vec<_> = (0..mesh.tris.len())
        .map(|cell| Cell {
            density_kg_m3: 2.0 + 10.0 * cell as f64,
            specific: [],
        })
        .collect();
    let u = [0.7, -0.2];
    let old_velocity = constant(&mesh, u);
    let velocity_boundary = |_| u;
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, face)| {
            (face.tris.1 == usize::MAX).then_some(Boundary::Velocity {
                edge,
                value: &velocity_boundary,
            })
        })
        .collect();
    let inlet: Vec<_> = mesh
        .edges
        .iter()
        .map(|face| {
            (face.tris.1 == usize::MAX).then_some(Cell {
                density_kg_m3: 321.0,
                specific: [],
            })
        })
        .collect();
    let exterior = vec![None; mesh.edges.len()];
    let mut expected_density = None;
    for reference in [0.5, 7.0, 1000.0] {
        let assembly = system
            .assemble_conservative_momentum_bdf1(
                History {
                    density: &density,
                    velocity: &old_velocity,
                    inflow_density: &inlet,
                    inflow_velocity: &exterior,
                },
                controls(reference),
                MaterialFluxNumerics::default(),
                &|_, _| 0.03,
                &|_| [0.0; 2],
                &boundaries,
                || false,
            )
            .unwrap();
        assert!(
            assembly
                .density()
                .iter()
                .zip(&density)
                .any(|(new, old)| new != old)
        );
        assert!(core::ptr::eq(assembly.mass_flux().mesh(), &mesh));
        let accepted = assembly.solve(|| false).unwrap();
        for (actual, expected) in accepted.field().x[..system.n_u].iter().zip(&old_velocity) {
            near(*actual, *expected, 3e-10);
        }
        for &pressure in accepted.pressure_pa() {
            near(pressure, 0.0, 3e-8);
        }
        near(
            accepted.mass_transport().relative_mass_balance_residual,
            0.0,
            3e-15,
        );
        let mass_per_depth = accepted.mass_transport().mass_after_kg / controls(reference).depth_m;
        near(
            accepted.kinetic_energy_j_per_m(),
            0.5 * mass_per_depth * (u[0] * u[0] + u[1] * u[1]),
            3e-9,
        );
        if let Some(expected) = &expected_density {
            assert_eq!(accepted.density(), expected);
        } else {
            expected_density = Some(accepted.density().to_vec());
        }
    }
}

#[test]
fn original_density_weights_physical_force_and_acceleration_under_natural_traction() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 1));
    let system = FluxSystem::new(&mesh);
    let rho = |point: [f64; 2]| if point[0] < 0.5 { 3.0 } else { 70.0 };
    let density: Vec<_> = mesh
        .centroids
        .iter()
        .map(|&point| Cell {
            density_kg_m3: rho(point),
            specific: [],
        })
        .collect();
    let old_velocity = vec![0.0; system.n_u];
    let inlet = vec![None; mesh.edges.len()];
    let exterior = vec![None; mesh.edges.len()];
    let traction = |_, _| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, face)| {
            (face.tris.1 == usize::MAX).then_some(Boundary::TotalTraction {
                edge,
                value: &traction,
            })
        })
        .collect();
    let acceleration = [1.7, -3.2];
    let force = |point| acceleration.map(|value| rho(point) * value);
    for reference in [1.0, 100.0] {
        let setting = controls(reference);
        let accepted = system
            .assemble_conservative_momentum_bdf1(
                History {
                    density: &density,
                    velocity: &old_velocity,
                    inflow_density: &inlet,
                    inflow_velocity: &exterior,
                },
                setting,
                MaterialFluxNumerics::default(),
                &|_, _| 0.2,
                &force,
                &boundaries,
                || false,
            )
            .unwrap()
            .solve(|| false)
            .unwrap();
        let expected = constant(&mesh, acceleration.map(|value| setting.timestep_s * value));
        for (actual, expected) in accepted.field().x[..system.n_u].iter().zip(expected) {
            near(*actual, expected, 2e-11);
        }
        assert_eq!(accepted.density(), density);
        let mass_per_depth = accepted.mass_transport().mass_after_kg / setting.depth_m;
        let actual_momentum = momentum(&system, accepted.density(), &accepted.field().x);
        for axis in 0..2 {
            near(
                actual_momentum[axis],
                mass_per_depth * setting.timestep_s * acceleration[axis],
                1e-10,
            );
        }
        for &pressure in accepted.pressure_pa() {
            near(pressure, 0.0, 2e-9);
        }
    }
}

#[test]
fn reversing_variable_density_flow_conserves_total_momentum_with_zero_boundary_exchange() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 0.8, 2, 2));
    let system = FluxSystem::new(&mesh);
    let density: Vec<_> = (0..mesh.tris.len())
        .map(|cell| Cell {
            density_kg_m3: 1.0 + 11.0 * cell as f64,
            specific: [],
        })
        .collect();
    let mut old_velocity = vec![0.0; system.n_u];
    for (edge, face) in mesh.edges.iter().enumerate() {
        if face.tris.1 != usize::MAX {
            old_velocity[2 * edge + 1] = 0.004 * (edge + 1) as f64;
        }
    }
    let inlet = vec![None; mesh.edges.len()];
    let exterior = vec![None; mesh.edges.len()];
    let traction = |_, _| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, face)| {
            (face.tris.1 == usize::MAX).then_some(Boundary::TotalTraction {
                edge,
                value: &traction,
            })
        })
        .collect();
    let before = momentum(&system, &density, &old_velocity);
    assert!(before[0].abs() + before[1].abs() > 1e-3);
    let accepted = system
        .assemble_conservative_momentum_bdf1(
            History {
                density: &density,
                velocity: &old_velocity,
                inflow_density: &inlet,
                inflow_velocity: &exterior,
            },
            controls(13.0),
            MaterialFluxNumerics::default(),
            &|_, _| 0.017,
            &|_| [0.0; 2],
            &boundaries,
            || false,
        )
        .unwrap()
        .solve(|| false)
        .unwrap();
    assert!(
        accepted
            .density()
            .iter()
            .zip(&density)
            .any(|(new, old)| new != old)
    );
    let after = momentum(&system, accepted.density(), &accepted.field().x);
    for axis in 0..2 {
        near(after[axis], before[axis], 2e-10);
    }
    near(
        accepted.mass_transport().mass_before_kg,
        accepted.mass_transport().mass_after_kg,
        3e-13,
    );
    assert_eq!(accepted.mass_transport().inlet_mass_kg, 0.0);
    assert_eq!(accepted.mass_transport().outlet_mass_kg, 0.0);
    for (edge, face) in mesh.edges.iter().enumerate() {
        if face.tris.1 == usize::MAX {
            assert_eq!(
                accepted.mass_flux().edge_directional_mass_rates(edge),
                Some([0.0; 2])
            );
        }
    }
}

#[test]
fn complete_step_refuses_late_cancellation_without_modifying_either_history() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let density = vec![
        Cell {
            density_kg_m3: 2.0,
            specific: [],
        },
        Cell {
            density_kg_m3: 9.0,
            specific: [],
        },
    ];
    let velocity = constant(&mesh, [0.3, -0.2]);
    let inlet: Vec<_> = mesh
        .edges
        .iter()
        .map(|face| (face.tris.1 == usize::MAX).then_some(density[0]))
        .collect();
    let exterior = vec![None; mesh.edges.len()];
    let prescribed = |_| [0.3, -0.2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, face)| {
            (face.tris.1 == usize::MAX).then_some(Boundary::Velocity {
                edge,
                value: &prescribed,
            })
        })
        .collect();
    let history = || History {
        density: &density,
        velocity: &velocity,
        inflow_density: &inlet,
        inflow_velocity: &exterior,
    };
    let mut calls = 0;
    let first = system
        .assemble_conservative_momentum_bdf1(
            history(),
            controls(1.0),
            MaterialFluxNumerics::default(),
            &|_, _| 0.1,
            &|_| [0.0; 2],
            &boundaries,
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        system.assemble_conservative_momentum_bdf1(
            history(),
            controls(1.0),
            MaterialFluxNumerics::default(),
            &|_, _| 0.1,
            &|_| [0.0; 2],
            &boundaries,
            || {
                seen += 1;
                seen == calls
            }
        ),
        Err(Error::Transport(DensityTransportError::Cancelled))
    ));
    assert_eq!(seen, calls);
    let mut solve_calls = 0;
    let accepted = first
        .solve(|| {
            solve_calls += 1;
            false
        })
        .unwrap();
    let second = system
        .assemble_conservative_momentum_bdf1(
            history(),
            controls(1.0),
            MaterialFluxNumerics::default(),
            &|_, _| 0.1,
            &|_| [0.0; 2],
            &boundaries,
            || false,
        )
        .unwrap();
    let mut solve_seen = 0;
    assert!(matches!(
        second.solve(|| {
            solve_seen += 1;
            solve_seen == solve_calls
        }),
        Err(Error::Transport(DensityTransportError::Cancelled))
    ));
    let retry = system
        .assemble_conservative_momentum_bdf1(
            history(),
            controls(1.0),
            MaterialFluxNumerics::default(),
            &|_, _| 0.1,
            &|_| [0.0; 2],
            &boundaries,
            || false,
        )
        .unwrap()
        .solve(|| false)
        .unwrap();
    assert_eq!(accepted.field().x, retry.field().x);
    assert_eq!(accepted.density(), retry.density());
    assert_eq!(density[0].density_kg_m3, 2.0);
    assert_eq!(density[1].density_kg_m3, 9.0);
    assert_eq!(velocity, constant(&mesh, [0.3, -0.2]));
}

#[path = "conservative_density_momentum/gravity.rs"]
mod gravity;
