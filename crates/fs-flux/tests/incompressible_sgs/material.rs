//! Original material rheology plus SGS at one physical nonlinear endpoint.
use super::*;
use fs_flux::material::{
    MaterialFluxError, MaterialFluxNumerics, MaterialFluxSystem, MaterialNavierStokesBoundary,
    MaterialSgsConfig, MaterialSgsError,
};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use std::cell::Cell;

struct Sample {
    state: FluidTransportState,
    rate: f64,
    cell: usize,
}
impl FluidTransport for Sample {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.state)
    }
}

// Manufactured shear thinning: the 0.6 factor has units seconds. This fixture
// checks coupling to an original source, not a calibrated constitutive claim.
fn molecular(rate: f64) -> f64 {
    0.1 + 0.5 / (1.0 + 0.6 * rate)
}
fn source(density: f64, rate: f64, cell: usize) -> Result<Sample, MaterialError> {
    Ok(Sample {
        state: FluidTransportState::try_new(density, 1000.0, 0.3, density * molecular(rate))?,
        rate,
        cell,
    })
}
fn controls() -> NavierStokesControls {
    NavierStokesControls {
        maximum_iterations: 120,
        relative_tolerance: 1e-10,
        relaxation: 0.8,
    }
}
fn numerics() -> MaterialFluxNumerics {
    MaterialFluxNumerics {
        linear_tolerance: 1e-12,
        ..MaterialFluxNumerics::default()
    }
}

#[test]
fn material_sgs_couples_actual_rheology_full_gradient_and_physical_heating_at_every_endpoint() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let density = 7.0;
    let reference = source(density, 0.0, usize::MAX).unwrap();
    let material = MaterialFluxSystem::new(&system, &reference, numerics()).unwrap();
    let widths = vec![0.7; mesh.tris.len()];
    let gradients = [[[0.3, 0.2], [0.0, -0.3]], [[0.3, 0.6], [-0.4, -0.3]]];
    for model in [
        SgsModel::Smagorinsky { cs: 0.5 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.07 },
    ] {
        let mut transfers = Vec::new();
        for (case, g) in gradients.into_iter().enumerate() {
            let target = |p: [f64; 2]| {
                [
                    0.4 + g[0][0] * p[0] + g[0][1] * p[1],
                    g[1][0] * p[0] + g[1][1] * p[1],
                ]
            };
            // Same symmetric strain, different rotation. Material rheology
            // must agree while the full-gradient models distinguish the cases.
            let rate = 0.4_f64.sqrt();
            let eddy = trace_free_planar_reference(model, 0.7, g);
            assert!(eddy > 0.0);
            let total = molecular(rate) + eddy;
            let traction = |p: [f64; 2], n: [f64; 2]| {
                let pressure = 5.0 + 2.0 * p[0] - 3.0 * p[1];
                let shear = total * (g[0][1] + g[1][0]);
                [
                    density * ((2.0 * total * g[0][0] - pressure) * n[0] + shear * n[1]),
                    density * (shear * n[0] + (2.0 * total * g[1][1] - pressure) * n[1]),
                ]
            };
            let slip_wall = |p: [f64; 2]| [target(p)[0] - 0.25 * total, 0.0];
            let boundaries: Vec<_> = mesh
                .edges
                .iter()
                .enumerate()
                .filter(|(edge, _)| mesh.is_boundary(*edge))
                .map(|(edge, face)| {
                    if face.normal[0].abs() > 0.5 {
                        MaterialNavierStokesBoundary::Open {
                            edge,
                            traction: &traction,
                            backflow: OpenBoundaryPolicy::NaturalTraction,
                        }
                    } else if case == 0 && face.normal[1] < 0.0 {
                        MaterialNavierStokesBoundary::NavierSlip {
                            edge,
                            value: &slip_wall,
                            tangential_resistance_pa_s_per_m: density * 0.8,
                        }
                    } else {
                        MaterialNavierStokesBoundary::Velocity {
                            edge,
                            value: &target,
                        }
                    }
                })
                .collect();
            let previous = moments(&mesh, &|p| target(p).map(|u| 0.8 * u));
            let older_equal = moments(&mesh, &|p| target(p).map(|u| 0.6 * u));
            let older_variable = moments(&mesh, &|p| target(p).map(|u| 0.5 * u));
            let saved = [
                previous.clone(),
                older_equal.clone(),
                older_variable.clone(),
            ];
            for scheme in 0..4 {
                // A(t)u_base with dA/dt=2 gives exactly the same derivative
                // for BDF1, equal BDF2 and dt=.1 / dt_old=.15 variable BDF2.
                let force = |p: [f64; 2]| {
                    let u = target(p);
                    let time = if scheme == 0 { 0.0 } else { 2.0 };
                    [
                        density * (time * u[0] + g[0][0] * u[0] + g[0][1] * u[1] + 2.0),
                        density * (time * u[1] + g[1][0] * u[0] + g[1][1] * u[1] - 3.0),
                    ]
                };
                let calls = Cell::new(0usize);
                let resolve = |point: fs_flux::material::MaterialFluxRheologyPoint,
                               _: &FluidTransportState| {
                    calls.set(calls.get() + 1);
                    source(density, point.shear_rate_per_s, point.cell)
                };
                let config = MaterialSgsConfig {
                    model,
                    filter_widths_m: &widths,
                };
                let step = match scheme {
                    0 => material.solve_navier_stokes_mixed_rheology_sgs(
                        config,
                        &force,
                        &boundaries,
                        &vec![0.0; system.n_u],
                        controls(),
                        resolve,
                    ),
                    1 => material.bdf1_navier_stokes_mixed_rheology_sgs_step(
                        config,
                        &force,
                        &boundaries,
                        &previous,
                        0.1,
                        controls(),
                        resolve,
                    ),
                    2 => material.bdf2_navier_stokes_mixed_rheology_sgs_step(
                        config,
                        &force,
                        &boundaries,
                        &previous,
                        Some(&older_equal),
                        0.1,
                        controls(),
                        resolve,
                    ),
                    _ => material.bdf2_navier_stokes_mixed_rheology_sgs_variable_step(
                        config,
                        &force,
                        &boundaries,
                        &previous,
                        Some((&older_variable, 0.15)),
                        0.1,
                        controls(),
                        resolve,
                    ),
                }
                .unwrap();
                let accepted_queries = calls.get();
                let flow = &step.flow.step;
                let (error, divergence) =
                    system.velocity_error(&flow.solution.discrete().x, &target);
                assert!(
                    error < 7e-9 && divergence < 7e-9,
                    "scheme {scheme}: {error}, {divergence}"
                );
                assert!(flow.relative_update <= controls().relative_tolerance);
                assert!(flow.nonlinear_residual <= controls().relative_tolerance);
                for (p, actual) in mesh.centroids.iter().zip(flow.solution.pressure_pa()) {
                    close(*actual, density * (5.0 + 2.0 * p[0] - 3.0 * p[1]), 7e-7);
                }
                for (cell, (sample, sgs)) in
                    step.flow.samples().iter().zip(&step.sgs.cells).enumerate()
                {
                    assert_eq!(sample.material().cell, cell);
                    close(sample.material().rate, rate, 5e-8);
                    close(
                        sample.point().shear_rate_per_s,
                        sgs.strain_rate_per_s,
                        1e-13,
                    );
                    close(
                        sample.transport().dynamic_viscosity_pa_s(),
                        density * molecular(rate),
                        5e-8,
                    );
                    close(sgs.eddy_viscosity_m2_s, eddy, 1e-8);
                    close(sgs.effective_viscosity_m2_s, total, 1e-8);
                }
                close(
                    flow.solution.viscous_dissipation().total_power_w_per_m(),
                    density * molecular(rate) * 0.4,
                    1e-7,
                );
                close(step.modeled_transfer_w_per_m, density * eddy * 0.4, 1e-7);
                assert_eq!(
                    calls.get(),
                    accepted_queries,
                    "reporting must retain original source owners"
                );
                assert!(accepted_queries >= mesh.tris.len());
                assert_eq!(previous, saved[0]);
                assert_eq!(older_equal, saved[1]);
                assert_eq!(older_variable, saved[2]);
                if scheme == 0 {
                    transfers.push(step.modeled_transfer_w_per_m);
                }
            }
        }
        if matches!(model, SgsModel::Smagorinsky { .. }) {
            close(transfers[0], transfers[1], 1e-8);
        } else {
            assert!((transfers[0] - transfers[1]).abs() > 1e-4);
        }
    }
}

#[test]
fn zero_material_sgs_preserves_native_variable_step_rheology_and_source_samples() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let density = 3.0;
    let reference = source(density, 0.0, usize::MAX).unwrap();
    let material = MaterialFluxSystem::new(&system, &reference, numerics()).unwrap();
    let widths = vec![0.5; mesh.tris.len()];
    let velocity = |p: [f64; 2]| [0.2 + 0.1 * p[0], -0.1 * p[1]];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| MaterialNavierStokesBoundary::Velocity {
            edge,
            value: &velocity,
        })
        .collect();
    let previous = moments(&mesh, &|p| velocity(p).map(|u| 0.8 * u));
    let older = moments(&mesh, &|p| velocity(p).map(|u| 0.5 * u));
    let force = |p: [f64; 2]| {
        let u = velocity(p);
        [density * 2.1 * u[0], density * 1.9 * u[1]]
    };
    let base = material
        .bdf2_navier_stokes_mixed_rheology_variable_step(
            &force,
            &boundaries,
            &previous,
            Some((&older, 0.15)),
            0.1,
            controls(),
            |point, _| source(density, point.shear_rate_per_s, point.cell),
        )
        .unwrap();
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let step = material
            .bdf2_navier_stokes_mixed_rheology_sgs_variable_step(
                MaterialSgsConfig {
                    model,
                    filter_widths_m: &widths,
                },
                &force,
                &boundaries,
                &previous,
                Some((&older, 0.15)),
                0.1,
                controls(),
                |point, _| source(density, point.shear_rate_per_s, point.cell),
            )
            .unwrap();
        assert_eq!(
            step.flow.step.solution.discrete().x,
            base.step.solution.discrete().x
        );
        assert_eq!(
            step.flow.step.solution.pressure_pa(),
            base.step.solution.pressure_pa()
        );
        assert_eq!(step.modeled_transfer_w_per_m, 0.0);
        for (sample, original) in step.flow.samples().iter().zip(base.samples()) {
            assert_eq!(sample.point(), original.point());
            assert_eq!(sample.transport(), original.transport());
        }
    }
}

#[test]
fn material_sgs_preserves_original_source_refusal_and_never_changes_history() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let density = 2.0;
    let reference = source(density, 0.0, usize::MAX).unwrap();
    let material = MaterialFluxSystem::new(&system, &reference, numerics()).unwrap();
    let widths = vec![0.5; mesh.tris.len()];
    let velocity = |p: [f64; 2]| [0.3 + 0.1 * p[1], 0.0];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| MaterialNavierStokesBoundary::Velocity {
            edge,
            value: &velocity,
        })
        .collect();
    let previous = moments(&mesh, &velocity);
    let saved = previous.clone();
    let config = MaterialSgsConfig {
        model: SgsModel::Smagorinsky { cs: 0.2 },
        filter_widths_m: &widths,
    };
    let calls = Cell::new(0);
    let refused = material.bdf1_navier_stokes_mixed_rheology_sgs_step(
        config,
        &|_| [0.0; 2],
        &boundaries,
        &previous,
        0.05,
        controls(),
        |point, _| {
            calls.set(calls.get() + 1);
            source(
                if point.cell + 1 == mesh.tris.len() {
                    2.0 * density
                } else {
                    density
                },
                point.shear_rate_per_s,
                point.cell,
            )
        },
    );
    assert!(matches!(
        refused,
        Err(MaterialSgsError::Material(MaterialFluxError::Invalid(_)))
    ));
    assert!(
        calls.get() > 1,
        "failure must occur after earlier source samples"
    );
    assert_eq!(previous, saved);
    let malformed = MaterialSgsConfig {
        model: SgsModel::Wale { cw: f64::NAN },
        ..config
    };
    calls.set(0);
    let refused = material.solve_navier_stokes_mixed_rheology_sgs(
        malformed,
        &|_| [0.0; 2],
        &boundaries,
        &previous,
        controls(),
        |point, _| {
            calls.set(calls.get() + 1);
            source(density, point.shear_rate_per_s, point.cell)
        },
    );
    assert!(matches!(
        refused,
        Err(MaterialSgsError::Subgrid(PlanarSgsError::InvalidInput(_)))
    ));
    assert_eq!(calls.get(), 0);
    assert_eq!(previous, saved);
}
