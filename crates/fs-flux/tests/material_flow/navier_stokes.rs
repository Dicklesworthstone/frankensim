//! Source density scales convection, physical traction and energy together.

use super::*;
use fs_flux::{MaterialNavierStokesBoundary, NavierStokesControls, OpenBoundaryPolicy};

#[path = "navier_stokes_rheology.rs"]
mod rheology;
#[path = "reservoir.rs"]
mod reservoir;

#[test]
fn nonlinear_mixed_flow_keeps_shared_gas_density_and_viscous_traction_in_si() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let spec = GasSpec::dry_air_ussa1976();
    let dt = 0.025;
    let target = |p: [f64; 2]| [0.2 * p[0] + 0.3 * p[1] + 0.4, -0.2 * p[1] + 0.1];
    let previous = vec![0.0; system.n_u];
    let mut energy_per_density = Vec::new();
    for source_pressure in [101_325.0, 202_650.0] {
        let gas = GasState::try_new(&spec, 293.15, source_pressure).unwrap();
        let owner =
            MaterialFluxSystem::new(&system, &gas, MaterialFluxNumerics::default()).unwrap();
        let density = owner.transport().density_kg_m3();
        let mu = owner.transport().dynamic_viscosity_pa_s();
        let force = |p: [f64; 2]| {
            let u = target(p);
            // grad(u)*u = [.04*x+.11, .04*y-.02].
            [
                density * (u[0] / dt + 0.04 * p[0] + 0.11) + 2.0,
                density * (u[1] / dt + 0.04 * p[1] - 0.02) - 3.0,
            ]
        };
        let traction = |p: [f64; 2], n: [f64; 2]| {
            let pressure = 7.0 + 2.0 * p[0] - 3.0 * p[1];
            [
                (0.4 * mu - pressure) * n[0] + 0.3 * mu * n[1],
                0.3 * mu * n[0] + (-0.4 * mu - pressure) * n[1],
            ]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(e, _)| mesh.is_boundary(*e))
            .map(|(edge, e)| {
                if e.normal[0].abs() > 0.5 {
                    MaterialNavierStokesBoundary::Open {
                        edge,
                        traction: &traction,
                        backflow: OpenBoundaryPolicy::NaturalTraction,
                    }
                } else {
                    MaterialNavierStokesBoundary::Velocity {
                        edge,
                        value: &target,
                    }
                }
            })
            .collect();
        let step = owner
            .bdf1_navier_stokes_mixed_step(
                &force,
                &boundaries,
                &previous,
                dt,
                NavierStokesControls {
                    relative_tolerance: 1e-10,
                    ..NavierStokesControls::default()
                },
            )
            .unwrap();
        let (error, divergence) = system.velocity_error(&step.solution.discrete().x, &target);
        assert!(error < 2e-9 && divergence < 2e-9, "{error}, {divergence}");
        for (actual, p) in step.solution.pressure_pa().iter().zip(&mesh.centroids) {
            near(*actual, 7.0 + 2.0 * p[0] - 3.0 * p[1], 3e-8);
        }
        assert_eq!(step.solution.transport(), owner.transport());
        assert_eq!(
            step.solution.pressure_reference(),
            StokesPressureReference::TractionDatum
        );
        assert!(step.iterations > 1 && step.nonlinear_residual < 1e-10);
        assert!(
            step.open_boundary_policies
                .iter()
                .all(|(_, p)| *p == OpenBoundaryPolicy::NaturalTraction)
        );
        energy_per_density.push(step.solution.kinetic_energy_j_per_m() / density);
    }
    // Independent integral on [0,1]² of |u|²/2 for the affine target.
    let x_energy = 0.2_f64.powi(2) / 3.0
        + 0.3_f64.powi(2) / 3.0
        + 0.4_f64.powi(2)
        + 0.2 * 0.3 / 2.0
        + 0.2 * 0.4
        + 0.3 * 0.4;
    let y_energy = 0.2_f64.powi(2) / 3.0 - 0.2 * 0.1 + 0.1_f64.powi(2);
    for value in energy_per_density {
        near(value, 0.5 * (x_energy + y_energy), 2e-9);
    }
    assert!(previous.iter().all(|value| *value == 0.0));
}

#[test]
fn steady_material_flow_has_live_convection_and_no_initial_history_load() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let target = |point: [f64; 2]| [0.4 + 0.3 * point[1], 0.1];
    // u.grad(u)=(0.03,0), gamma=0.3 and integral(|u|^2)/2=0.16.
    for density in [2.0, 7.0] {
        let mu = 0.5 * density;
        let source = transport(density, mu);
        let owner =
            MaterialFluxSystem::new(&system, &source, MaterialFluxNumerics::default()).unwrap();
        let force = |_: [f64; 2]| [0.03 * density + 2.0, -3.0];
        let traction = |point: [f64; 2], normal: [f64; 2]| {
            let pressure = 7.0 + 2.0 * point[0] - 3.0 * point[1];
            [
                -pressure * normal[0] + 0.3 * mu * normal[1],
                0.3 * mu * normal[0] - pressure * normal[1],
            ]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, geometry)| {
                if geometry.normal[0].abs() > 0.5 {
                    MaterialNavierStokesBoundary::Open {
                        edge,
                        traction: &traction,
                        backflow: OpenBoundaryPolicy::NaturalTraction,
                    }
                } else {
                    MaterialNavierStokesBoundary::Velocity {
                        edge,
                        value: &target,
                    }
                }
            })
            .collect();
        for initial_speed in [0.0, 0.2] {
            let initial: Vec<_> = mesh
                .edges
                .iter()
                .flat_map(|edge| [initial_speed * edge.normal[0], 0.0])
                .collect();
            let snapshot = initial.clone();
            let result = owner
                .solve_navier_stokes_mixed(
                    &force,
                    &boundaries,
                    &initial,
                    NavierStokesControls {
                        relative_tolerance: 1e-10,
                        ..NavierStokesControls::default()
                    },
                )
                .unwrap();
            let (error, divergence) = system.velocity_error(&result.solution.discrete().x, &target);
            assert!(error < 3e-9 && divergence < 3e-9, "{error}, {divergence}");
            assert!(result.nonlinear_residual <= 1e-10);
            assert_eq!(result.solution.transport(), &source);
            assert_eq!(
                result.solution.pressure_reference(),
                StokesPressureReference::TractionDatum
            );
            for (pressure, point) in result.solution.pressure_pa().iter().zip(&mesh.centroids) {
                near(*pressure, 7.0 + 2.0 * point[0] - 3.0 * point[1], 4e-8);
            }
            near(
                result.solution.kinetic_energy_j_per_m(),
                0.16 * density,
                3e-8,
            );
            near(
                result.solution.viscous_dissipation().total_power_w_per_m(),
                0.09 * mu,
                3e-8,
            );
            assert_eq!(initial, snapshot);
        }
    }
}

#[test]
fn variable_step_material_flow_retains_physical_pressure_inertia_and_heat() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let dt = 0.1;
    let time = 0.3;
    let amplitude = |t: f64| 1.0 + 0.2 * t + 0.3 * t * t;
    let scale = amplitude(time);
    let base = |point: [f64; 2]| [0.4 + 0.3 * point[1], 0.1];
    let target = |point| base(point).map(|u| scale * u);
    let base_moments: Vec<_> = mesh
        .edges
        .iter()
        .flat_map(|edge| {
            let u = base(edge.mid);
            let dy = mesh.verts[edge.verts.1][1] - mesh.verts[edge.verts.0][1];
            [
                u[0] * edge.normal[0] + u[1] * edge.normal[1],
                0.3 * dy * edge.normal[0] / 12.0,
            ]
        })
        .collect();
    let previous: Vec<_> = base_moments
        .iter()
        .map(|u| amplitude(time - dt) * u)
        .collect();
    for previous_dt in [0.05, 0.2] {
        let older: Vec<_> = base_moments
            .iter()
            .map(|u| amplitude(time - dt - previous_dt) * u)
            .collect();
        let snapshots = [previous.clone(), older.clone()];
        for density in [2.0, 7.0] {
            let mu = 0.4;
            let source = transport(density, mu);
            let owner =
                MaterialFluxSystem::new(&system, &source, MaterialFluxNumerics::default()).unwrap();
            let force = |point| {
                let v = base(point);
                [
                    density * ((0.2 + 0.6 * time) * v[0] + 0.03 * scale * scale) + 2.0,
                    density * (0.2 + 0.6 * time) * v[1] - 3.0,
                ]
            };
            let traction = |point: [f64; 2], normal: [f64; 2]| {
                let pressure = 7.0 + 2.0 * point[0] - 3.0 * point[1];
                [
                    -pressure * normal[0] + 0.3 * scale * mu * normal[1],
                    0.3 * scale * mu * normal[0] - pressure * normal[1],
                ]
            };
            let boundaries: Vec<_> = mesh
                .edges
                .iter()
                .enumerate()
                .filter(|(edge, _)| mesh.is_boundary(*edge))
                .map(|(edge, geometry)| {
                    if geometry.normal[0].abs() > 0.5 {
                        MaterialNavierStokesBoundary::Open {
                            edge,
                            traction: &traction,
                            backflow: OpenBoundaryPolicy::NaturalTraction,
                        }
                    } else {
                        MaterialNavierStokesBoundary::Velocity {
                            edge,
                            value: &target,
                        }
                    }
                })
                .collect();
            let result = owner
                .bdf2_navier_stokes_mixed_variable_step(
                    &force,
                    &boundaries,
                    &previous,
                    Some((&older, previous_dt)),
                    dt,
                    NavierStokesControls {
                        relative_tolerance: 1e-10,
                        ..NavierStokesControls::default()
                    },
                )
                .unwrap();
            let solution = &result.solution;
            let (error, divergence) = system.velocity_error(&solution.discrete().x, &target);
            assert!(error < 3e-9 && divergence < 3e-9, "{error}, {divergence}");
            for (pressure, point) in solution.pressure_pa().iter().zip(&mesh.centroids) {
                near(*pressure, 7.0 + 2.0 * point[0] - 3.0 * point[1], 4e-8);
            }
            near(
                solution.kinetic_energy_j_per_m(),
                0.16 * density * scale * scale,
                4e-8,
            );
            near(
                solution.viscous_dissipation().total_power_w_per_m(),
                mu * 0.09 * scale * scale,
                4e-8,
            );
            assert_eq!(solution.transport(), &source);
            assert_eq!(previous, snapshots[0]);
            assert_eq!(older, snapshots[1]);
        }
    }
}
