#![cfg(feature = "continuum")]
//! G1/G3: density transport, native BDM momentum and physical boundary exchange
//! use one compatible mass flux, including sign-changing faces.

use fs_flux::bdm::edge_gauss_pub;
use fs_flux::ns::boundary::conservative::{
    DensityMomentumControls as Controls, DensityMomentumError as Error,
};
use fs_flux::ns::transport::density::{DensityTransportCell as Cell, DensityTransportError};
use fs_flux::{
    FluxParams, FluxSystem, MaterialFluxNumerics, MaterialNavierStokesBoundary as Boundary,
    NavierStokesBoundary, OpenBoundaryPolicy, StokesBoundary, TriMesh,
};
use fs_solid::Mesh2;

#[path = "density_momentum/subcycle.rs"]
mod subcycle;

#[path = "density_momentum/ssprk2.rs"]
mod ssprk2;

#[cfg(feature = "experimental-sgs")]
#[path = "density_momentum/sgs.rs"]
mod sgs;

fn controls(reference: f64, depth: f64) -> Controls {
    Controls {
        reference_density_kg_m3: reference,
        depth_m: depth,
        maximum_courant: 0.8,
        relative_volume_divergence_tolerance: 1e-12,
        maximum_transport_storage_bytes: 1 << 22,
    }
}

fn numerics() -> MaterialFluxNumerics {
    MaterialFluxNumerics {
        penalty: 30.0,
        linear_tolerance: 1e-11,
        max_linear_iterations: 20_000,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn moments(system: &FluxSystem<'_>, field: &dyn Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    let mesh = system.mesh;
    let mut result = vec![0.0; system.n_u];
    for (index, edge) in mesh.edges.iter().enumerate() {
        let a = mesh.verts[edge.verts.0];
        let b = mesh.verts[edge.verts.1];
        for (point, weight) in edge_gauss_pub(a, b) {
            let u = field(point);
            let normal = u[0] * edge.normal[0] + u[1] * edge.normal[1];
            let s = ((point[0] - a[0]) * (b[0] - a[0]) + (point[1] - a[1]) * (b[1] - a[1]))
                / (edge.len * edge.len)
                - 0.5;
            result[2 * index] += weight * normal / edge.len;
            result[2 * index + 1] += weight * normal * s / edge.len;
        }
    }
    result
}

#[test]
fn transported_density_preserves_native_free_stream_and_physical_boundary_momentum() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 1));
    let system = FluxSystem::new(&mesh);
    let u = [0.4, -0.2];
    let uniform = |_: [f64; 2]| u;
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -7.0 * n);
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| Boundary::OpenReservoir {
            edge,
            traction: &traction,
            velocity: &uniform,
        })
        .collect();
    let selected: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, _)| mesh.is_boundary(edge).then_some(edge))
        .collect();
    let initial_velocity = moments(&system, &uniform);
    for (low, high) in [(2.0, 8.0), (2e-8, 8e8)] {
        let initial: Vec<Cell<1>> = mesh
            .centroids
            .iter()
            .map(|point| Cell {
                density_kg_m3: if point[0] < 0.5 { low } else { high },
                specific: [123.0],
            })
            .collect();
        let inflow: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .map(|(edge, _)| {
                mesh.is_boundary(edge).then_some(Cell {
                    density_kg_m3: high,
                    specific: [123.0],
                })
            })
            .collect();
        let mut normalized = None;
        for reference in [high, 0.125 * high] {
            let step = system
                .density_momentum_bdf1_step(
                    numerics(),
                    &|_, _, _| 0.05 * high,
                    &|_, _, _| [0.0; 2],
                    &boundaries,
                    &initial_velocity,
                    &initial_velocity,
                    &initial,
                    &inflow,
                    0.04,
                    controls(reference, 2.0),
                    || false,
                )
                .unwrap();
            let (error, divergence) = system.velocity_error(&step.flow.field.x, &uniform);
            assert!(
                error < 2e-10 && divergence < 2e-10,
                "{error:e}, {divergence:e}"
            );
            for &pressure in &step.pressure_pa {
                // Absolute pressure resolution follows the inertial scale
                // rho*U/dt; the extreme-contrast case must still retain its datum.
                near(pressure, 7.0, 2e-9 + 4096.0 * f64::EPSILON * high / 0.04);
            }
            assert!(
                step.cells
                    .iter()
                    .zip(&initial)
                    .any(|(a, b)| a.density_kg_m3 != b.density_kg_m3)
            );
            assert!(step.cells.iter().all(|cell| cell.specific == [123.0]));
            let scale = step
                .transport
                .mass_before_kg
                .max(step.transport.mass_after_kg)
                .max(1.0);
            near(step.transport.mass_balance_residual_kg, 0.0, 2e-14 * scale);
            for axis in 0..2 {
                near(
                    step.observations.momentum_after_kg_m_s[axis],
                    step.transport.mass_after_kg * u[axis],
                    3e-10 * scale,
                );
            }
            let out = step
                .flow
                .reactions
                .transport(&selected, [0.0; 2], [0.0; 2])
                .unwrap();
            near(out.outward_volume_flux_m2_s, 0.0, 2e-15);
            for axis in 0..2 {
                near(
                    step.observations.momentum_after_kg_m_s[axis]
                        - step.observations.momentum_before_kg_m_s[axis]
                        + 0.04 * reference * 2.0 * out.outward_momentum.force[axis],
                    0.0,
                    3e-10 * scale,
                );
            }
            near(
                step.observations.kinetic_energy_after_j
                    - step.observations.kinetic_energy_before_j
                    + 0.04 * reference * 2.0 * out.outward_kinetic_energy_flux,
                0.0,
                3e-10 * scale,
            );
            assert!(step.observations.temporal_dissipation_j < 1e-18 * scale);
            if let Some(previous) = normalized {
                for (a, b) in step.flow.field.x[..system.n_u].iter().zip(previous) {
                    near(*a, b, 2e-10);
                }
            }
            normalized = Some(step.flow.field.x[..system.n_u].to_vec());
        }
    }
}

#[test]
fn closed_unequal_density_counterflow_closes_native_energy_and_wall_momentum() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let mut previous = vec![0.0; system.n_u];
    previous[2 * internal + 1] = 0.2 / (1.5 * mesh.edges[internal].len);
    let cells = [
        Cell::<0> {
            density_kg_m3: 2.0,
            specific: [],
        },
        Cell {
            density_kg_m3: 8.0,
            specific: [],
        },
    ];
    let zero = |_: [f64; 2]| [0.0; 2];
    let selected: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, _)| mesh.is_boundary(edge).then_some(edge))
        .collect();
    let boundaries: Vec<_> = selected
        .iter()
        .map(|&edge| Boundary::Velocity { edge, value: &zero })
        .collect();
    let native: Vec<_> = selected
        .iter()
        .map(|&edge| StokesBoundary::Velocity { edge, value: &zero })
        .collect();
    let reference = 4.0;
    let depth = 1.7;
    let dt = 0.05;
    let mu = 0.03;
    let viscosity = system
        .assemble_stokes_mixed_variable_viscosity(
            FluxParams {
                nu: mu / reference,
                sigma: numerics().penalty,
                ..FluxParams::default()
            },
            &|_, _| mu / reference,
            &zero,
            &native,
        )
        .unwrap();
    let step = system
        .density_momentum_bdf1_step(
            numerics(),
            &|_, _, _| mu,
            &|_, _, _| [0.0; 2],
            &boundaries,
            &previous,
            &previous,
            &cells,
            &vec![None; mesh.edges.len()],
            dt,
            controls(reference, depth),
            || false,
        )
        .unwrap();
    near(
        step.transport.mass_before_kg,
        step.transport.mass_after_kg,
        3e-14,
    );
    assert!(step.cells[0].density_kg_m3 > cells[0].density_kg_m3);
    assert!(step.cells[1].density_kg_m3 < cells[1].density_kg_m3);
    let mut velocity = step.flow.field.x.clone();
    velocity[system.n_u..].fill(0.0);
    let mut action = vec![0.0; system.n];
    viscosity.matrix().spmv(&velocity, &mut action);
    let viscous_work: f64 = velocity
        .iter()
        .zip(action)
        .map(|(u, au)| u * au)
        .sum::<f64>()
        * dt
        * reference
        * depth;
    let obs = step.observations;
    assert!(obs.temporal_dissipation_j > 0.0 && obs.interior_upwind_dissipation_j > 0.0);
    assert!(viscous_work > 0.0);
    near(
        obs.kinetic_energy_after_j - obs.kinetic_energy_before_j
            + obs.temporal_dissipation_j
            + obs.interior_upwind_dissipation_j
            + viscous_work,
        0.0,
        2e-12,
    );
    let loads = step
        .flow
        .reactions
        .resultant(&selected, [0.0; 2], &[Some(0.0)])
        .unwrap();
    for axis in 0..2 {
        near(
            obs.momentum_after_kg_m_s[axis] - obs.momentum_before_kg_m_s[axis]
                + dt * reference * depth * loads.stress.force[axis],
            0.0,
            2e-12,
        );
    }
}

#[test]
fn physical_force_uses_actual_cell_density_and_retains_traction_pressure() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 1));
    let system = FluxSystem::new(&mesh);
    let old = vec![0.0; system.n_u];
    let cells: Vec<Cell<0>> = mesh
        .centroids
        .iter()
        .map(|p| Cell {
            density_kg_m3: if p[0] < 0.5 { 2.0 } else { 9.0 },
            specific: [],
        })
        .collect();
    let traction = |_: [f64; 2], normal: [f64; 2]| normal.map(|n| -11.0 * n);
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| Boundary::Open {
            edge,
            traction: &traction,
            backflow: OpenBoundaryPolicy::NaturalTraction,
        })
        .collect();
    let acceleration = [0.7, -0.3];
    let dt = 0.2;
    let calls = std::cell::Cell::new(0);
    let step = system
        .density_momentum_bdf1_step(
            numerics(),
            &|_, _, state| 0.03 * state.density_kg_m3,
            &|_, _, state| {
                calls.set(calls.get() + 1);
                acceleration.map(|a| a * state.density_kg_m3)
            },
            &boundaries,
            &old,
            &old,
            &cells,
            &vec![None; mesh.edges.len()],
            dt,
            controls(0.3, 1.0),
            || false,
        )
        .unwrap();
    assert_eq!(calls.get(), 6 * cells.len());
    let (error, divergence) =
        system.velocity_error(&step.flow.field.x, &|_| acceleration.map(|a| dt * a));
    assert!(error < 5e-11 && divergence < 5e-11);
    for p in &step.pressure_pa {
        near(*p, 11.0, 3e-10);
    }
    assert_eq!(step.cells, cells);
    for axis in 0..2 {
        near(
            step.observations.momentum_after_kg_m_s[axis],
            step.transport.mass_after_kg * dt * acceleration[axis],
            2e-10,
        );
    }
}

/// G1: an independent closed-form solution of the spatially discrete coupled
/// system. There is one free velocity moment m on the internal diagonal. Its
/// affine basis has integral |phi|²=12 in EACH triangle, squared tangential
/// jump 144, and directional volume rate (3*sqrt(2)/2)*m. Thus total density
/// S stays fixed, m'=-lambda*m-beta*m² and D'=-6*sqrt(2)*m*D, where
/// The symmetric-stress SIP form is 96*mu*(4*sigma-3), using the normal
/// cell height in its penalty. Therefore lambda=96*mu*(4*sigma-3)/(12*S),
/// beta=9*sqrt(2), and D=rho_left-rho_right.
#[test]
fn coupled_density_and_momentum_have_first_order_temporal_convergence() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| Boundary::Velocity { edge, value: &zero })
        .collect();
    let inflow = vec![None; mesh.edges.len()];
    let mu = 0.02;
    let initial_m = 0.1;
    let end = 0.2;
    let lambda = 96.0 * mu * (4.0 * numerics().penalty - 3.0) / 120.0;
    let beta = 9.0 * 2.0f64.sqrt();
    let denominator = 1.0 + beta * initial_m / lambda * (-f64::exp_m1(-lambda * end));
    let exact_m = initial_m * f64::exp(-lambda * end) / denominator;
    let exact_difference = -6.0 * denominator.powf(-2.0 / 3.0);
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let dt = end / f64::from(steps);
        let mut previous = vec![0.0; system.n_u];
        previous[2 * internal + 1] = initial_m;
        let mut cells = vec![
            Cell::<0> {
                density_kg_m3: 2.0,
                specific: [],
            },
            Cell {
                density_kg_m3: 8.0,
                specific: [],
            },
        ];
        for _ in 0..steps {
            let step = system
                .density_momentum_bdf1_step(
                    numerics(),
                    &|_, _, _| mu,
                    &|_, _, _| [0.0; 2],
                    &boundaries,
                    &previous,
                    &previous,
                    &cells,
                    &inflow,
                    dt,
                    controls(5.0, 1.0),
                    || false,
                )
                .unwrap();
            previous = step.flow.field.x[..system.n_u].to_vec();
            cells = step.cells;
        }
        near(cells[0].density_kg_m3 + cells[1].density_kg_m3, 10.0, 8e-14);
        errors.push([
            (previous[2 * internal + 1] - exact_m).abs(),
            (cells[0].density_kg_m3 - cells[1].density_kg_m3 - exact_difference).abs(),
        ]);
    }
    for pair in errors.windows(2) {
        for component in 0..2 {
            let ratio = pair[0][component] / pair[1][component];
            assert!(
                (1.8..2.2).contains(&ratio),
                "component={component}, ratio={ratio:e}, errors={errors:?}"
            );
        }
    }
}

#[test]
fn uniform_density_recovers_original_native_oseen_operator() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let mut old = vec![0.0; system.n_u];
    old[2 * internal + 1] = 0.02;
    let cells = vec![
        Cell::<0> {
            density_kg_m3: 3.0,
            specific: []
        };
        mesh.tris.len()
    ];
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundary_edges: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter_map(|(edge, _)| mesh.is_boundary(edge).then_some(edge))
        .collect();
    let boundaries: Vec<_> = boundary_edges
        .iter()
        .map(|&edge| Boundary::Velocity { edge, value: &zero })
        .collect();
    let original: Vec<_> = boundary_edges
        .iter()
        .map(|&edge| NavierStokesBoundary::Velocity { edge, value: &zero })
        .collect();
    let native = system
        .assemble_navier_stokes_mixed_bdf1_variable_viscosity(
            FluxParams {
                nu: 0.05,
                sigma: numerics().penalty,
                ..FluxParams::default()
            },
            &|_, _| 0.05,
            &zero,
            &original,
            &old,
            &old,
            0.1,
        )
        .unwrap();
    let transported = system
        .assemble_density_momentum_bdf1(
            numerics(),
            &|_, _, _| 0.15,
            &|_, _, _| [0.0; 2],
            &boundaries,
            &old,
            &old,
            &cells,
            &vec![None; mesh.edges.len()],
            0.1,
            controls(3.0, 1.0),
            || false,
        )
        .unwrap();
    let a = native.matrix().to_dense();
    let b = transported.momentum().matrix().to_dense();
    for (x, y) in a.into_iter().zip(b) {
        near(x, y, 3e-12);
    }
    for (&x, &y) in native.rhs().iter().zip(transported.momentum().rhs()) {
        near(x, y, 3e-13);
    }
}

#[test]
fn cancellation_cfl_and_invalid_source_never_mutate_caller_state() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let internal = mesh
        .edges
        .iter()
        .position(|edge| edge.tris.1 != usize::MAX)
        .unwrap();
    let mut previous = vec![0.0; system.n_u];
    previous[2 * internal + 1] = 0.2;
    let cells = vec![
        Cell::<0> {
            density_kg_m3: 1.0,
            specific: []
        };
        2
    ];
    let retained_cells = cells.clone();
    let retained_velocity = previous.clone();
    let zero = |_: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| Boundary::Velocity { edge, value: &zero })
        .collect();
    let inflow = vec![None; mesh.edges.len()];
    let calls = std::cell::Cell::new(0);
    let force = |_, _, _: &Cell<0>| {
        calls.set(calls.get() + 1);
        [0.0; 2]
    };
    let refused = system.assemble_density_momentum_bdf1(
        numerics(),
        &|_, _, _| 0.1,
        &force,
        &boundaries,
        &previous,
        &previous,
        &cells,
        &inflow,
        100.0,
        controls(1.0, 1.0),
        || false,
    );
    assert!(matches!(
        refused,
        Err(Error::Transport(DensityTransportError::Courant { .. }))
    ));
    assert_eq!(calls.get(), 0);
    for stop in [1, 4, 15] {
        let mut polls = 0;
        let refused = system.assemble_density_momentum_bdf1(
            numerics(),
            &|_, _, _| 0.1,
            &force,
            &boundaries,
            &previous,
            &previous,
            &cells,
            &inflow,
            0.01,
            controls(1.0, 1.0),
            || {
                polls += 1;
                polls >= stop
            },
        );
        assert!(matches!(
            refused,
            Err(Error::Cancelled | Error::Transport(DensityTransportError::Cancelled))
        ));
    }
    let refused = system.assemble_density_momentum_bdf1(
        numerics(),
        &|_, _, _| -0.1,
        &force,
        &boundaries,
        &previous,
        &previous,
        &cells,
        &inflow,
        0.01,
        controls(1.0, 1.0),
        || false,
    );
    assert!(matches!(refused, Err(Error::Flow(_))));
    let assembly = system
        .assemble_density_momentum_bdf1(
            numerics(),
            &|_, _, _| 0.1,
            &force,
            &boundaries,
            &previous,
            &previous,
            &cells,
            &inflow,
            0.01,
            controls(1.0, 1.0),
            || false,
        )
        .unwrap();
    assert!(matches!(assembly.solve(|| true), Err(Error::Cancelled)));
    let quiescent = vec![0.0; system.n_u];
    let enormous_history = moments(&system, &|_| [1e200, 0.0]);
    let viscosity_calls = std::cell::Cell::new(0);
    let refused = system.assemble_density_momentum_bdf1(
        numerics(),
        &|_, _, _| {
            viscosity_calls.set(viscosity_calls.get() + 1);
            0.1
        },
        &force,
        &boundaries,
        &quiescent,
        &enormous_history,
        &cells,
        &inflow,
        1e-200,
        controls(1.0, 1.0),
        || false,
    );
    assert!(matches!(
        refused,
        Err(Error::Invalid(
            "unrepresentable physical momentum history or load"
        ))
    ));
    assert_eq!(viscosity_calls.get(), 0);
    assert_eq!(cells, retained_cells);
    assert_eq!(previous, retained_velocity);
}
