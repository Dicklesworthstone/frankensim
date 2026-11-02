//! Pressure-driven inertia and endpoint boundary data on the original FEEC law.

use super::{separate_domains, velocity_boundaries};
use fs_flux::bdm::edge_gauss_pub;
use fs_flux::{
    FluxParams, FluxSystem, StokesBoundary, StokesBoundaryError, StokesPressureReference, TriMesh,
};
use fs_solid::Mesh2;
use std::cell::Cell;

fn traction_boundaries<'a>(
    mesh: &TriMesh,
    value: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
) -> Vec<StokesBoundary<'a>> {
    mesh.edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::TotalTraction { edge, value })
        .collect()
}

/// With no velocity walls, a pressure gradient gives u' = -grad(p).
/// Viscosity cannot damp uniform translation; ending the pulse leaves coasting.
#[test]
fn pressure_pulse_uses_physical_inertia_without_a_steady_velocity_pin() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let mut previous = vec![0.0; system.n_u];
    let mut speed = 0.0;
    for (step, (dt, acceleration)) in [(0.07, 2.0), (0.13, 2.0), (0.11, 0.0), (0.05, -1.0)]
        .into_iter()
        .enumerate()
    {
        let datum = 4.0 + step as f64;
        let traction = |point: [f64; 2], normal: [f64; 2]| {
            let pressure = datum - acceleration * point[0];
            [-pressure * normal[0], -pressure * normal[1]]
        };
        let boundaries = traction_boundaries(&mesh, &traction);
        assert!(matches!(
            system.assemble_stokes_mixed(params, &|_| [0.0; 2], &boundaries),
            Err(StokesBoundaryError::NoVelocityBoundary)
        ));
        let snapshot = previous.clone();
        let assembly = system
            .assemble_stokes_mixed_bdf1(params, &|_| [0.0; 2], &boundaries, &previous, dt)
            .unwrap();
        let solved = assembly.solve().unwrap();
        assert_eq!(previous, snapshot);
        let old_speed = speed;
        speed += dt * acceleration;
        let (error, divergence) = system.velocity_error(&solved.field.x, &|_| [speed, 0.0]);
        assert!(error < 2e-11 && divergence < 2e-11, "{error}, {divergence}");
        assert_eq!(
            solved.pressure_reference,
            StokesPressureReference::TractionDatum
        );
        for (actual, point) in solved.pressure.iter().zip(&mesh.centroids) {
            assert!((actual - (datum - acceleration * point[0])).abs() < 2e-10);
        }
        // Backward-Euler work = kinetic change + its temporal dissipation.
        let kinetic_change =
            system.kinetic_energy(&solved.field.x) - system.kinetic_energy(&previous);
        let temporal_dissipation = 0.5 * (speed - old_speed).powi(2);
        let pressure_work = dt * acceleration * speed;
        assert!((kinetic_change + temporal_dissipation - pressure_work).abs() < 2e-12);
        if step == 0 {
            let sparse = system.solve_sparse(assembly.matrix(), assembly.rhs(), params);
            assert!(sparse.converged(params.tol));
            for (a, b) in sparse.x.iter().zip(&solved.field.x) {
                assert!((a - b).abs() < 3e-8);
            }
        }
        previous.copy_from_slice(&solved.field.x[..system.n_u]);
    }
}

/// p(x,t)=3+5t-t*x accelerates the open chamber with u'=t. The independent
/// continuous solution u=t²/2 exposes first-order temporal error exactly.
#[test]
fn pressure_ramp_is_first_order_and_keeps_component_local_datums() {
    let mesh = separate_domains();
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let final_time = 0.2;
    let mut errors = Vec::new();
    for steps in [2, 4, 8] {
        let dt = final_time / steps as f64;
        let mut previous = vec![0.0; system.n_u];
        for step in 1..=steps {
            let time = step as f64 * dt;
            let datum = 3.0 + 5.0 * time;
            let zero = |_: [f64; 2]| [0.0; 2];
            let traction = |point: [f64; 2], normal: [f64; 2]| {
                let pressure = datum - time * (point[0] - 3.0);
                [-pressure * normal[0], -pressure * normal[1]]
            };
            let boundaries: Vec<_> = mesh
                .edges
                .iter()
                .enumerate()
                .filter(|(edge, _)| mesh.is_boundary(*edge))
                .map(|(edge, geometry)| {
                    if geometry.mid[0] < 2.0 {
                        StokesBoundary::Velocity { edge, value: &zero }
                    } else {
                        StokesBoundary::TotalTraction {
                            edge,
                            value: &traction,
                        }
                    }
                })
                .collect();
            let solution = system
                .bdf1_stokes_mixed_step(
                    params,
                    &|point| {
                        if point[0] < 2.0 {
                            [1.3, -0.4]
                        } else {
                            [0.0; 2]
                        }
                    },
                    &boundaries,
                    &previous,
                    dt,
                )
                .unwrap();
            let discrete_speed = 0.5 * time * (time + dt);
            let (error, divergence) = system.velocity_error(&solution.field.x, &|point| {
                if point[0] < 2.0 {
                    [0.0; 2]
                } else {
                    [discrete_speed, 0.0]
                }
            });
            assert!(error < 3e-11 && divergence < 3e-11);
            assert_eq!(
                solution.pressure_reference,
                StokesPressureReference::Componentwise
            );
            for (cell, point) in mesh.centroids.iter().enumerate() {
                let (expected, reference) = if point[0] < 2.0 {
                    (
                        1.3 * (point[0] - 0.5) - 0.4 * (point[1] - 0.5),
                        StokesPressureReference::MeanZero,
                    )
                } else {
                    (
                        datum - time * (point[0] - 3.0),
                        StokesPressureReference::TractionDatum,
                    )
                };
                assert_eq!(solution.pressure_reference_per_cell[cell], reference);
                assert!((solution.pressure[cell] - expected).abs() < 2e-10);
            }
            previous.copy_from_slice(&solution.field.x[..system.n_u]);
        }
        let measured = system.velocity_at(&previous, [3.4, 0.4])[0];
        let error = measured - 0.5 * final_time * final_time;
        assert!((error - 0.5 * final_time * dt).abs() < 3e-11);
        errors.push(error);
    }
    assert!((errors[0] / errors[1] - 2.0).abs() < 1e-7);
    assert!((errors[1] / errors[2] - 2.0).abs() < 1e-7);
}

#[test]
fn mixed_bdf1_reuses_the_original_mass_and_symmetric_stress_operator() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let viscosity = |_: usize, point: [f64; 2]| 0.5 + point[1];
    let force = |point: [f64; 2]| [point[1], -point[0]];
    let velocity = |point: [f64; 2]| [0.2 * point[1], 0.0];
    let boundaries = velocity_boundaries(&mesh, &velocity);
    let previous: Vec<_> = (0..system.n_u).map(|i| 0.001 * i as f64).collect();
    let dt: f64 = 0.03;
    let mixed = system
        .assemble_stokes_mixed_bdf1_variable_viscosity(
            params,
            &viscosity,
            &force,
            &boundaries,
            &previous,
            dt,
        )
        .unwrap();
    let (matrix, rhs) = system.assemble_viscous_transient(
        params,
        &viscosity,
        &force,
        &velocity,
        None,
        dt.recip(),
        Some(&previous),
    );
    assert_eq!(mixed.matrix().to_dense(), matrix.to_dense());
    assert_eq!(mixed.rhs(), rhs.as_slice());
    assert_eq!(
        mixed.pressure_reference(),
        StokesPressureReference::MeanZero
    );
}

/// An exact affine shear with dnu/dy != 0 checks implicit inertia and the
/// viscous force together, including the total stress on both pressure ends.
#[test]
fn transient_variable_viscosity_matches_an_independent_affine_solution() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let viscosity = |_: usize, point: [f64; 2]| 0.7 + 1.2 * point[1];
    let dt = 0.2;
    let previous_velocity = |point: [f64; 2]| [0.3 + 0.4 * point[1], 0.0];
    // u(t,y) = 0.3 + 1.1*t + (0.4 + 0.6*t)*y; p=8-2x.
    let slope = 0.4 + 0.6 * dt;
    let velocity = |point: [f64; 2]| [0.3 + 1.1 * dt + slope * point[1], 0.0];
    let force = |point: [f64; 2]| [1.1 + 0.6 * point[1] - 1.2 * slope - 2.0, 0.0];
    let traction = |point: [f64; 2], normal: [f64; 2]| {
        let pressure = 8.0 - 2.0 * point[0];
        let shear = viscosity(0, point) * slope;
        [
            shear * normal[1] - pressure * normal[0],
            shear * normal[0] - pressure * normal[1],
        ]
    };
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, geometry)| {
            if geometry.normal[1].abs() > 0.5 {
                StokesBoundary::Velocity {
                    edge,
                    value: &velocity,
                }
            } else {
                StokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let mut previous = vec![0.0; system.n_u];
    for (edge, geometry) in mesh.edges.iter().enumerate() {
        let a = mesh.verts[geometry.verts.0];
        let b = mesh.verts[geometry.verts.1];
        for (point, weight) in edge_gauss_pub(a, b) {
            let u = previous_velocity(point);
            let un = u[0] * geometry.normal[0] + u[1] * geometry.normal[1];
            let s = ((point[0] - a[0]) * (b[0] - a[0]) + (point[1] - a[1]) * (b[1] - a[1]))
                / geometry.len.powi(2)
                - 0.5;
            previous[2 * edge] += weight * un / geometry.len;
            previous[2 * edge + 1] += weight * un * s / geometry.len;
        }
    }
    let solution = system
        .assemble_stokes_mixed_bdf1_variable_viscosity(
            params,
            &viscosity,
            &force,
            &boundaries,
            &previous,
            dt,
        )
        .unwrap()
        .solve()
        .unwrap();
    let (error, divergence) = system.velocity_error(&solution.field.x, &velocity);
    assert!(error < 2e-11 && divergence < 2e-11, "{error}, {divergence}");
    for (pressure, point) in solution.pressure.iter().zip(&mesh.centroids) {
        assert!((pressure - (8.0 - 2.0 * point[0])).abs() < 2e-10);
    }
}

#[test]
fn invalid_history_and_timestep_refuse_before_boundary_or_load_queries() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let params = FluxParams::default();
    let calls = Cell::new(0);
    let query = |_: [f64; 2]| {
        calls.set(calls.get() + 1);
        [0.0; 2]
    };
    let boundaries = velocity_boundaries(&mesh, &query);
    let previous = vec![0.0; system.n_u];
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(matches!(
            system.assemble_stokes_mixed_bdf1(params, &query, &boundaries, &previous, dt,),
            Err(StokesBoundaryError::InvalidInput(_))
        ));
    }
    for invalid in [vec![0.0; system.n_u - 1], vec![f64::NAN; system.n_u]] {
        assert!(matches!(
            system.assemble_stokes_mixed_bdf1(params, &query, &boundaries, &invalid, 0.1,),
            Err(StokesBoundaryError::InvalidInput(_))
        ));
    }
    assert_eq!(calls.get(), 0);
    assert!(previous.iter().all(|value| *value == 0.0));

    let leaking_velocity = |point: [f64; 2]| point;
    let boundaries = velocity_boundaries(&mesh, &leaking_velocity);
    assert!(matches!(
        system.assemble_stokes_mixed_bdf1(params, &query, &boundaries, &previous, 0.1,),
        Err(StokesBoundaryError::IncompatibleFlux { .. })
    ));
    assert_eq!(calls.get(), 0);
}
