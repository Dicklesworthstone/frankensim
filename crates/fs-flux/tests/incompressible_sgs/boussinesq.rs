//! Thermal feedback, reference-pressure balance and original endpoint oracles.

use super::*;
use fs_flux::ns::sgs::{
    PlanarBoussinesqConfig, PlanarSgsBoussinesqError, PlanarSgsThermalControls,
    PlanarSgsThermalInput,
};
use fs_flux::ns::transport::scalar::{ScalarBoundary, ScalarSolveControls};

fn config() -> PlanarBoussinesqConfig {
    PlanarBoussinesqConfig {
        reference_temperature_k: 305.0,
        thermal_expansion_per_k: 0.002,
        gravity_m_s2: [0.0, -9.81],
        maximum_relative_density_departure: 0.05,
    }
}
fn controls() -> PlanarSgsThermalControls {
    PlanarSgsThermalControls {
        flow: NavierStokesControls {
            maximum_iterations: 100,
            relative_tolerance: 1e-9,
            relaxation: 0.8,
        },
        thermal_solve: ScalarSolveControls {
            relative_tolerance: 1e-12,
            maximum_iterations: 1000,
        },
        ..PlanarSgsThermalControls::default()
    }
}
fn params() -> FluxParams {
    FluxParams {
        nu: f64::NAN,
        tol: 1e-13,
        ..FluxParams::default()
    }
}
fn values(mesh: &TriMesh, f: &dyn Fn([f64; 2]) -> f64) -> Vec<f64> {
    mesh.tris
        .iter()
        .flat_map(|tri| tri.map(|v| f(mesh.verts[v])))
        .collect()
}
fn thermal_boundaries<'a>(
    mesh: &TriMesh,
    value: &'a dyn Fn([f64; 2]) -> f64,
) -> Vec<ScalarBoundary<'a>> {
    mesh.edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, _)| ScalarBoundary::Value { edge, value })
        .collect()
}
fn sample(mesh: &TriMesh, temperature: &[f64], p: [f64; 2]) -> f64 {
    let cell = mesh.containing_cell(p).unwrap();
    let [a, b, c] = mesh.tris[cell].map(|v| mesh.verts[v]);
    let area2 = 2.0 * mesh.areas[cell];
    let phi = [
        ((b[0] - p[0]) * (c[1] - p[1]) - (b[1] - p[1]) * (c[0] - p[0])) / area2,
        ((c[0] - p[0]) * (a[1] - p[1]) - (c[1] - p[1]) * (a[0] - p[0])) / area2,
        ((a[0] - p[0]) * (b[1] - p[1]) - (a[1] - p[1]) * (b[0] - p[0])) / area2,
    ];
    (0..3).map(|i| phi[i] * temperature[3 * cell + i]).sum()
}

#[test]
fn all_closures_evolve_temperature_buoyancy_and_momentum_at_one_endpoint() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3));
    let system = FluxSystem::new(&mesh);
    let wall = |_: [f64; 2]| [0.0; 2];
    let boundary_temperature = |p: [f64; 2]| 300.0 + 10.0 * p[0];
    let source = |_: usize, _: [f64; 2]| 0.0;
    let force = |_: [f64; 2]| [0.0; 2];
    let flow_bc = velocity_boundaries(&mesh, &wall);
    let thermal_bc = thermal_boundaries(&mesh, &boundary_temperature);
    let old_u = vec![0.0; system.n_u];
    let old_t = values(&mesh, &boundary_temperature);
    let snapshot = old_t.clone();
    let alpha = vec![0.02; mesh.tris.len()];
    let prandtl = vec![0.7; mesh.tris.len()];
    let input = PlanarSgsThermalInput {
        previous_temperature_k: &old_t,
        molecular_diffusivity_m2_s: &alpha,
        turbulent_prandtl: &prandtl,
        boundaries: &thermal_bc,
        source_k_s: &source,
    };
    for closure in [
        SgsModel::Smagorinsky { cs: 0.3 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.1 },
    ] {
        let model = PlanarSgs::with_model(
            &system,
            closure,
            &vec![0.4; mesh.tris.len()],
            &vec![0.03; mesh.tris.len()],
        )
        .unwrap();
        let dt = 0.2;
        let result = model
            .bdf1_boussinesq_step(
                params(),
                config(),
                &force,
                &flow_bc,
                &old_u,
                input,
                dt,
                controls(),
                || false,
            )
            .unwrap();
        let u = &result.flow.solution.field.x[..system.n_u];
        let warm_side_up: f64 = mesh
            .centroids
            .iter()
            .enumerate()
            .map(|(cell, &p)| mesh.areas[cell] * (p[0] - 0.5) * system.velocity_at(u, p)[1])
            .sum();
        assert!(warm_side_up > 1e-5, "hot side must rise: {warm_side_up:e}");
        assert!(result.buoyancy_power_m4_s3 > 0.0);
        assert!(
            result
                .temperature
                .values
                .iter()
                .zip(&old_t)
                .any(|(a, b)| (a - b).abs() > 1e-5)
        );
        assert!(
            result
                .sgs
                .cells
                .iter()
                .any(|c| c.eddy_viscosity_m2_s > 1e-10)
        );
        assert!(result.flow.nonlinear_residual <= controls().flow.relative_tolerance);
        assert!(
            result.temperature.relative_residual <= controls().thermal_solve.relative_tolerance
        );
        close(result.temperature.balance_residual, 0.0, 2e-8);
        // An independent native thermal re-solve at the accepted velocity
        // must reproduce the returned endpoint, not an earlier Picard state.
        let native = model
            .prepare_scalar_transport(
                u,
                &alpha,
                &prandtl,
                &thermal_bc,
                controls().transport,
                || false,
            )
            .unwrap();
        let replay = native
            .transport
            .bdf1_step(&old_t, dt, &source, controls().thermal_solve, || false)
            .unwrap();
        assert_eq!(replay.values, result.temperature.values);
        assert_eq!(native.sgs, result.sgs);

        let final_force = |p| {
            [
                0.0,
                0.002 * 9.81 * (sample(&mesh, &result.temperature.values, p) - 305.0),
            ]
        };
        let fixed_temperature = model
            .bdf1_step(
                params(),
                &final_force,
                &flow_bc,
                &old_u,
                dt,
                controls().flow,
            )
            .unwrap();
        for (a, b) in u
            .iter()
            .zip(&fixed_temperature.step.solution.field.x[..system.n_u])
        {
            close(*a, *b, 3e-9);
        }
        let lagged_force = |p: [f64; 2]| [0.0, 0.002 * 9.81 * (boundary_temperature(p) - 305.0)];
        let lagged = model
            .bdf1_step(
                params(),
                &lagged_force,
                &flow_bc,
                &old_u,
                dt,
                controls().flow,
            )
            .unwrap();
        assert!(
            u.iter()
                .zip(&lagged.step.solution.field.x)
                .any(|(a, b)| (a - b).abs() > 1e-8)
        );

        let passive = model
            .bdf1_boussinesq_step(
                params(),
                PlanarBoussinesqConfig {
                    thermal_expansion_per_k: 0.0,
                    ..config()
                },
                &force,
                &flow_bc,
                &old_u,
                input,
                dt,
                controls(),
                || false,
            )
            .unwrap();
        assert!(
            passive.flow.solution.field.x[..system.n_u]
                .iter()
                .all(|&v| v == 0.0)
        );
        close(passive.buoyancy_power_m4_s3, 0.0, 0.0);
    }
    assert_eq!(old_t, snapshot);
    assert!(old_u.iter().all(|&v| v == 0.0));
}

#[test]
fn uniform_heating_uses_endpoint_buoyancy_and_closes_the_original_thermal_budget() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let velocity = |_: [f64; 2]| [0.0, 0.1];
    let old_u = moments(&mesh, &velocity);
    let flow_bc = velocity_boundaries(&mesh, &velocity);
    let donor = |_: [f64; 2]| 300.2;
    let zero = |_: [f64; 2]| 0.0;
    let thermal_bc: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, _)| ScalarBoundary::Flux {
            edge,
            outward: &zero,
            inflow: Some(&donor),
        })
        .collect();
    let old_t = vec![300.0; 3 * mesh.tris.len()];
    let alpha = vec![0.02; mesh.tris.len()];
    let prandtl = vec![0.7; mesh.tris.len()];
    let source = |_: usize, _: [f64; 2]| 2.0;
    let config = PlanarBoussinesqConfig {
        reference_temperature_k: 300.0,
        thermal_expansion_per_k: 0.01,
        gravity_m_s2: [0.0, -10.0],
        maximum_relative_density_departure: 0.05,
    };
    let model = PlanarSgs::new(
        &system,
        0.3,
        &vec![0.25; mesh.tris.len()],
        &vec![0.05; mesh.tris.len()],
    )
    .unwrap();
    let result = model
        .bdf1_boussinesq_step(
            params(),
            config,
            &|_| [0.0, -0.02],
            &flow_bc,
            &old_u,
            PlanarSgsThermalInput {
                previous_temperature_k: &old_t,
                molecular_diffusivity_m2_s: &alpha,
                turbulent_prandtl: &prandtl,
                boundaries: &thermal_bc,
                source_k_s: &source,
            },
            0.1,
            controls(),
            || false,
        )
        .unwrap();
    for value in &result.temperature.values {
        close(*value, 300.2, 1e-9);
    }
    for (a, b) in result.flow.solution.field.x[..system.n_u]
        .iter()
        .zip(&old_u)
    {
        close(*a, *b, 1e-10);
    }
    close(result.temperature.storage_rate, 2.0, 1e-8);
    close(result.temperature.source_rate, 2.0, 1e-12);
    close(result.temperature.balance_residual, 0.0, 1e-8);
    close(result.buoyancy_power_m4_s3, 0.002, 1e-11);
    close(result.maximum_relative_density_departure, 0.002, 2e-11);
}

#[test]
fn signed_gravity_retains_reference_reduced_hydrostatic_pressure_and_no_spurious_heat() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let velocity = |_: [f64; 2]| [0.2, 0.0];
    let old_u = moments(&mesh, &velocity);
    let flow_bc = velocity_boundaries(&mesh, &velocity);
    let donor = |_: [f64; 2]| 310.0;
    let thermal_bc = thermal_boundaries(&mesh, &donor);
    let old_t = vec![310.0; 3 * mesh.tris.len()];
    let alpha = vec![0.02; mesh.tris.len()];
    let prandtl = vec![0.7; mesh.tris.len()];
    let source = |_: usize, _: [f64; 2]| 0.0;
    let model = PlanarSgs::new(
        &system,
        0.3,
        &vec![0.25; mesh.tris.len()],
        &vec![0.05; mesh.tris.len()],
    )
    .unwrap();
    for gravity in [-10.0, 10.0] {
        let result = model
            .bdf1_boussinesq_step(
                params(),
                PlanarBoussinesqConfig {
                    reference_temperature_k: 300.0,
                    thermal_expansion_per_k: 0.001,
                    gravity_m_s2: [0.0, gravity],
                    maximum_relative_density_departure: 0.05,
                },
                &|_| [0.0; 2],
                &flow_bc,
                &old_u,
                PlanarSgsThermalInput {
                    previous_temperature_k: &old_t,
                    molecular_diffusivity_m2_s: &alpha,
                    turbulent_prandtl: &prandtl,
                    boundaries: &thermal_bc,
                    source_k_s: &source,
                },
                0.1,
                controls(),
                || false,
            )
            .unwrap();
        for (a, b) in result.flow.solution.field.x[..system.n_u]
            .iter()
            .zip(&old_u)
        {
            close(*a, *b, 1e-10);
        }
        for (cell, p) in mesh.centroids.iter().enumerate() {
            close(
                result.flow.solution.field.x[system.n_u + cell],
                // The original closed owner pins its FIRST pressure cell.
                -0.01 * gravity * (p[1] - mesh.centroids[0][1]),
                1e-9,
            );
        }
        for value in &result.temperature.values {
            close(*value, 310.0, 1e-9);
        }
        close(result.temperature.source_rate, 0.0, 0.0);
        close(result.temperature.balance_residual, 0.0, 1e-8);
        close(result.buoyancy_power_m4_s3, 0.0, 1e-11);
    }
}

#[test]
fn failed_coupling_or_late_cancellation_never_publishes_temperature_or_velocity() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let velocity = |_: [f64; 2]| [0.1, 0.0];
    let old_u = moments(&mesh, &velocity);
    let flow_bc = velocity_boundaries(&mesh, &velocity);
    let donor = |p: [f64; 2]| 300.0 + 10.0 * p[0];
    let thermal_bc = thermal_boundaries(&mesh, &donor);
    let old_t = values(&mesh, &donor);
    let snapshot_t = old_t.clone();
    let snapshot_u = old_u.clone();
    let alpha = vec![0.02; mesh.tris.len()];
    let prandtl = vec![0.7; mesh.tris.len()];
    let source = |_: usize, _: [f64; 2]| 0.0;
    let input = PlanarSgsThermalInput {
        previous_temperature_k: &old_t,
        molecular_diffusivity_m2_s: &alpha,
        turbulent_prandtl: &prandtl,
        boundaries: &thermal_bc,
        source_k_s: &source,
    };
    let model = PlanarSgs::new(
        &system,
        0.3,
        &vec![0.25; mesh.tris.len()],
        &vec![0.05; mesh.tris.len()],
    )
    .unwrap();
    let mut polls = 0;
    let expected = model
        .bdf1_boussinesq_step(
            params(),
            config(),
            &|_| [0.0; 2],
            &flow_bc,
            &old_u,
            input,
            0.1,
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut observed = 0;
    assert!(matches!(
        model.bdf1_boussinesq_step(
            params(),
            config(),
            &|_| [0.0; 2],
            &flow_bc,
            &old_u,
            input,
            0.1,
            controls(),
            || {
                observed += 1;
                observed >= polls
            }
        ),
        Err(PlanarSgsBoussinesqError::Cancelled)
    ));
    let exhausted = PlanarSgsThermalControls {
        flow: NavierStokesControls {
            maximum_iterations: 1,
            ..controls().flow
        },
        ..controls()
    };
    assert!(matches!(
        model.bdf1_boussinesq_step(
            params(),
            config(),
            &|_| [0.0; 2],
            &flow_bc,
            &old_u,
            input,
            0.1,
            exhausted,
            || false
        ),
        Err(PlanarSgsBoussinesqError::Flow(_))
    ));
    let bad_source = |_: usize, _: [f64; 2]| f64::NAN;
    assert!(matches!(
        model.bdf1_boussinesq_step(
            params(),
            config(),
            &|_| [0.0; 2],
            &flow_bc,
            &old_u,
            PlanarSgsThermalInput {
                source_k_s: &bad_source,
                ..input
            },
            0.1,
            controls(),
            || false
        ),
        Err(PlanarSgsBoussinesqError::Thermal(_))
    ));
    assert!(matches!(
        model.bdf1_boussinesq_step(
            params(),
            PlanarBoussinesqConfig {
                maximum_relative_density_departure: 0.11,
                ..config()
            },
            &|_| [0.0; 2],
            &flow_bc,
            &old_u,
            input,
            0.1,
            controls(),
            || false
        ),
        Err(PlanarSgsBoussinesqError::Invalid(_))
    ));
    let invalid_t = vec![400.0; old_t.len()];
    assert!(matches!(
        model.bdf1_boussinesq_step(
            params(),
            config(),
            &|_| [0.0; 2],
            &flow_bc,
            &old_u,
            PlanarSgsThermalInput {
                previous_temperature_k: &invalid_t,
                ..input
            },
            0.1,
            controls(),
            || false
        ),
        Err(PlanarSgsBoussinesqError::Temperature { .. })
    ));
    // Both Gauss samples meet the density-departure limit, but the affine
    // trace extrapolated to the edge's left endpoint does not. Admission must
    // cover the actual entire P1 boundary donor, not just quadrature samples.
    let edge = mesh
        .edges
        .iter()
        .enumerate()
        .find(|(e, geometry)| {
            mesh.is_boundary(*e)
                && mesh.verts[geometry.verts.0][1] == 0.0
                && mesh.verts[geometry.verts.1][1] == 0.0
                && mesh.verts[geometry.verts.0][0].min(mesh.verts[geometry.verts.1][0]) == 0.0
        })
        .unwrap()
        .0;
    let trace = |p: [f64; 2]| 275.0 + 100.0 * p[0];
    let mut bad_boundary = thermal_bc.clone();
    for boundary in &mut bad_boundary {
        if matches!(boundary, ScalarBoundary::Value { edge: e, .. } if *e == edge) {
            *boundary = ScalarBoundary::Value {
                edge,
                value: &trace,
            };
        }
    }
    assert!(matches!(
        model.bdf1_boussinesq_step(
            params(),
            config(),
            &|_| [0.0; 2],
            &flow_bc,
            &old_u,
            PlanarSgsThermalInput {
                boundaries: &bad_boundary,
                ..input
            },
            0.1,
            controls(),
            || false,
        ),
        Err(PlanarSgsBoussinesqError::BoundaryTemperature { .. })
    ));
    let retry = model
        .bdf1_boussinesq_step(
            params(),
            config(),
            &|_| [0.0; 2],
            &flow_bc,
            &old_u,
            input,
            0.1,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(expected.temperature.values, retry.temperature.values);
    assert_eq!(expected.flow.solution.field.x, retry.flow.solution.field.x);
    assert_eq!(old_t, snapshot_t);
    assert_eq!(old_u, snapshot_u);
}
