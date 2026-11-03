//! Analytic wall momentum and energy, through the original BDM1/P0 owners.

use fs_flux::ale::translating::TranslatingAleFlow;
use fs_flux::{
    BoundaryReactionOwner, FluxParams, FluxSystem, MaterialFluxNumerics, MaterialFluxSystem,
    MaterialNavierStokesBoundary, MaterialStokesBoundary, NavierStokesBoundary,
    NavierStokesControls, OpenBoundaryPolicy, StokesBoundary, StokesBoundaryError, TriMesh,
};
use fs_material::fluid_transport::FluidTransportState;
use fs_solid::Mesh2;
use std::cell::Cell;

fn mesh() -> TriMesh {
    TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3))
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "{actual} versus {expected}"
    );
}
fn zero(_: [f64; 2]) -> [f64; 2] {
    [0.0; 2]
}
fn controls() -> NavierStokesControls {
    NavierStokesControls {
        relative_tolerance: 1e-10,
        ..NavierStokesControls::default()
    }
}
// Independent exact edge integrals of an affine vector field.
fn affine_history(mesh: &TriMesh, field: impl Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    let mut values = vec![0.0; 2 * mesh.edges.len()];
    for (edge, geometry) in mesh.edges.iter().enumerate() {
        let a = field(mesh.verts[geometry.verts.0]);
        let b = field(mesh.verts[geometry.verts.1]);
        for axis in 0..2 {
            values[2 * edge] += 0.5 * (a[axis] + b[axis]) * geometry.normal[axis];
            values[2 * edge + 1] += (b[axis] - a[axis]) * geometry.normal[axis] / 12.0;
        }
    }
    values
}
fn selected(mesh: &TriMesh, side: f64) -> Vec<usize> {
    mesh.edges
        .iter()
        .enumerate()
        .filter(|(e, g)| mesh.is_boundary(*e) && g.normal[1] * side > 0.5)
        .map(|(e, _)| e)
        .collect()
}

/// Exact Couette solution: tau = nu*(U1-U0)/(H+nu/beta0+nu/beta1).
/// Its independent wall work equals bulk shear loss plus both slip losses.
#[test]
fn finite_and_perfect_slip_recover_couette_stress_torque_and_physical_loss() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let params = FluxParams {
        nu: 0.7,
        ..FluxParams::default()
    };
    let bottom = selected(&mesh, -1.0);
    let top = selected(&mesh, 1.0);
    let all: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .collect();
    let beta1 = 2.1;
    let (u0, u1, pressure) = (0.2, 1.1, 3.0);
    for beta0 in [0.0, 1.3] {
        let tau = if beta0 == 0.0 {
            0.0
        } else {
            params.nu * (u1 - u0) / (1.0 + params.nu / beta0 + params.nu / beta1)
        };
        let intercept = if beta0 == 0.0 { u1 } else { u0 + tau / beta0 };
        let exact = |x: [f64; 2]| [intercept + tau / params.nu * x[1], 0.0];
        let calls = Cell::new(0);
        let wall = |x: [f64; 2]| {
            calls.set(calls.get() + 1);
            [if x[1] < 0.5 { u0 } else { u1 }, 0.0]
        };
        let traction = |_: [f64; 2], n: [f64; 2]| {
            [-pressure * n[0] + tau * n[1], tau * n[0] - pressure * n[1]]
        };
        let boundaries: Vec<_> = all
            .iter()
            .map(|&edge| {
                let n = mesh.edges[edge].normal;
                if n[1].abs() > 0.5 {
                    StokesBoundary::NavierSlip {
                        edge,
                        value: &wall,
                        tangential_resistance_m_s: if n[1] < 0.0 { beta0 } else { beta1 },
                    }
                } else {
                    StokesBoundary::TotalTraction {
                        edge,
                        value: &traction,
                    }
                }
            })
            .collect();
        let result = system
            .assemble_stokes_mixed(params, &zero, &boundaries)
            .unwrap()
            .solve()
            .unwrap();
        let (error, div) = system.velocity_error(&result.field.x, &exact);
        assert!(error < 2e-9 && div < 2e-9);
        assert_eq!(calls.get(), 2 * (bottom.len() + top.len()));
        let origin = [0.2, -0.4];
        let low = result
            .reactions
            .resultant(&bottom, origin, &[None])
            .unwrap();
        let high = result.reactions.resultant(&top, origin, &[None]).unwrap();
        near(low.stress.force[0], tau, 3e-8);
        near(low.stress.force[1], -pressure, 3e-8);
        near(low.stress.moment, -0.3 * pressure - 0.4 * tau, 3e-8);
        near(high.stress.force[0], -tau, 3e-8);
        near(high.stress.force[1], pressure, 3e-8);
        let expected_loss = if beta0 == 0.0 { 0.0 } else { tau * tau / beta0 };
        near(low.slip_dissipation, expected_loss, 2e-9);
        near(high.slip_dissipation, tau * tau / beta1, 2e-9);
        let supplied_work = -(low.stress.force[0] * u0 + high.stress.force[0] * u1);
        near(
            supplied_work,
            tau * tau / params.nu + low.slip_dissipation + high.slip_dissipation,
            3e-8,
        );
        let total = result.reactions.resultant(&all, origin, &[None]).unwrap();
        near(total.stress.force[0], 0.0, 3e-8);
        near(total.stress.force[1], 0.0, 3e-8);
        near(total.stress.moment, 0.0, 3e-8);
        for item in result.reactions.edges() {
            if bottom.contains(&item.edge) || top.contains(&item.edge) {
                assert_eq!(item.owner, BoundaryReactionOwner::NavierSlip);
            }
        }
        assert_eq!(
            calls.get(),
            2 * (bottom.len() + top.len()),
            "reporting must not resample walls"
        );
        let nonlinear: Vec<_> = all
            .iter()
            .map(|&edge| {
                let n = mesh.edges[edge].normal;
                if n[1].abs() > 0.5 {
                    NavierStokesBoundary::NavierSlip {
                        edge,
                        value: &wall,
                        tangential_resistance_m_s: if n[1] < 0.0 { beta0 } else { beta1 },
                    }
                } else {
                    NavierStokesBoundary::Open {
                        edge,
                        traction: &traction,
                        backflow: OpenBoundaryPolicy::NaturalTraction,
                    }
                }
            })
            .collect();
        let step = system
            .bdf1_navier_stokes_mixed_step(
                params,
                &zero,
                &nonlinear,
                &affine_history(&mesh, exact),
                0.05,
                controls(),
            )
            .unwrap();
        let (error, div) = system.velocity_error(&step.solution.field.x, &exact);
        assert!(error < 2e-8 && div < 2e-8);
        let dynamic = step
            .solution
            .reactions
            .resultant(&bottom, origin, &[None])
            .unwrap();
        near(dynamic.stress.force[0], tau, 5e-8);
        near(dynamic.slip_dissipation, expected_loss, 2e-8);
        assert_eq!(
            calls.get(),
            4 * (bottom.len() + top.len()),
            "Picard retains admitted wall samples"
        );
        // Nonlinear viscosity must feed the same wall/reaction operator.
        let rate = tau / params.nu;
        let viscosity = |_: usize, _: [f64; 2], gamma: f64| {
            params.nu * (1.0 + gamma * gamma) / (1.0 + rate * rate)
        };
        let rheology = system
            .bdf1_navier_stokes_mixed_rheology_step(
                params,
                &viscosity,
                &zero,
                &nonlinear,
                &affine_history(&mesh, exact),
                0.05,
                controls(),
            )
            .unwrap();
        let (error, div) = system.velocity_error(&rheology.solution.field.x, &exact);
        assert!(error < 3e-8 && div < 3e-8);
        let load = rheology
            .solution
            .reactions
            .resultant(&bottom, origin, &[None])
            .unwrap();
        near(load.stress.force[0], tau, 6e-8);
        near(load.slip_dissipation, expected_loss, 3e-8);
        assert_eq!(calls.get(), 6 * (bottom.len() + top.len()));
    }
}

/// Density scales inertia, beta and viscosity together; a common translation
/// and acceleration preserve slip traction and dissipated power.
#[test]
fn physical_resistance_survives_source_density_and_accelerating_ale_frames() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let (mu, beta, wall_speed) = (0.4, 1.6, 0.3);
    let tau = mu * wall_speed / (1.0 + 2.0 * mu / beta);
    let exact = |x: [f64; 2]| [tau / beta + tau / mu * x[1], 0.0];
    let bottom = selected(&mesh, -1.0);
    let top = selected(&mesh, 1.0);
    let wall = |x: [f64; 2]| [if x[1] < 0.5 { 0.0 } else { wall_speed }, 0.0];
    let traction = |_: [f64; 2], n: [f64; 2]| [tau * n[1], tau * n[0]];
    for rho in [2.0, 7.0] {
        let material = FluidTransportState::try_new(rho, 1000.0, 0.5, mu).unwrap();
        let flow =
            MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(e, _)| mesh.is_boundary(*e))
            .map(|(edge, g)| {
                if g.normal[1].abs() > 0.5 {
                    MaterialStokesBoundary::NavierSlip {
                        edge,
                        value: &wall,
                        tangential_resistance_pa_s_per_m: beta,
                    }
                } else {
                    MaterialStokesBoundary::TotalTraction {
                        edge,
                        value: &traction,
                    }
                }
            })
            .collect();
        let solved = flow
            .bdf1_stokes_mixed_step(&zero, &boundaries, &affine_history(&mesh, exact), 0.05)
            .unwrap();
        let load = solved.boundary_loads(&bottom, [0.0; 2], &[None]).unwrap();
        near(load.stress.force_n_per_m[0], tau, 5e-8);
        near(load.slip_dissipation_w_per_m, tau * tau / beta, 2e-8);
        for boost in [[0.0; 2], [0.7, -0.4]] {
            let initial = affine_history(&mesh, |x| {
                let u = exact(x);
                [u[0] + boost[0], u[1] + boost[1]]
            });
            let mut owner = TranslatingAleFlow::new(
                &system,
                &material,
                MaterialFluxNumerics::default(),
                boost,
                &initial,
            )
            .unwrap();
            let dt = 0.05;
            let acceleration = [0.6, -0.4];
            let speed = [
                boost[0] + dt * acceleration[0],
                boost[1] + dt * acceleration[1],
            ];
            let low_wall = |_: [f64; 2]| speed;
            let high_wall = |_: [f64; 2]| [speed[0] + wall_speed, speed[1]];
            let moving: Vec<_> = mesh
                .edges
                .iter()
                .enumerate()
                .filter(|(e, _)| mesh.is_boundary(*e))
                .map(|(edge, g)| {
                    if g.normal[1].abs() > 0.5 {
                        let value: &dyn Fn([f64; 2]) -> [f64; 2] = if g.normal[1] < 0.0 {
                            &low_wall
                        } else {
                            &high_wall
                        };
                        MaterialNavierStokesBoundary::NavierSlip {
                            edge,
                            value,
                            tangential_resistance_pa_s_per_m: beta,
                        }
                    } else {
                        MaterialNavierStokesBoundary::Open {
                            edge,
                            traction: &traction,
                            backflow: OpenBoundaryPolicy::NaturalTraction,
                        }
                    }
                })
                .collect();
            let force = |_: [f64; 2]| acceleration.map(|a| rho * a);
            let step = owner
                .step(dt, speed, 1e-12, &force, &moving, controls(), || false)
                .unwrap();
            let (error, div) = system.velocity_error(&step.world_velocity_dofs, &|x| {
                let u = exact(x);
                [u[0] + speed[0], u[1] + speed[1]]
            });
            assert!(error < 2e-8 && div < 2e-8);
            let origin = step.translation_m;
            let low = step.boundary_loads_world(&bottom, origin, &[None]).unwrap();
            let high = step.boundary_loads_world(&top, origin, &[None]).unwrap();
            near(low.stress.force_n_per_m[0], tau, 5e-8);
            near(high.stress.force_n_per_m[0], -tau, 5e-8);
            near(low.slip_dissipation_w_per_m, tau * tau / beta, 2e-8);
            let transport = step.boundary_transport_world(&bottom, origin).unwrap();
            near(transport.outward_mass_flux_kg_per_m_s, 0.0, 2e-9);
        }
    }
}

#[test]
fn invalid_resistance_and_incompatible_normal_flux_refuse_before_solving() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let calls = Cell::new(0);
    let counted = |_: [f64; 2]| {
        calls.set(calls.get() + 1);
        [0.0; 2]
    };
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(e, _)| mesh.is_boundary(*e))
            .map(|(edge, _)| StokesBoundary::NavierSlip {
                edge,
                value: &counted,
                tangential_resistance_m_s: bad,
            })
            .collect();
        assert!(matches!(
            system.assemble_stokes_mixed(FluxParams::default(), &zero, &boundaries),
            Err(StokesBoundaryError::InvalidInput(_))
        ));
        assert_eq!(calls.get(), 0);
    }
    let material = FluidTransportState::try_new(2.0, 1000.0, 0.5, 0.4).unwrap();
    let flow =
        MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
    // A nonzero physical coefficient may not silently become perfect slip
    // after density normalization, including a tiny negative input.
    for beta in [f64::from_bits(1), -f64::from_bits(1)] {
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(e, _)| mesh.is_boundary(*e))
            .map(|(edge, _)| MaterialStokesBoundary::NavierSlip {
                edge,
                value: &counted,
                tangential_resistance_pa_s_per_m: beta,
            })
            .collect();
        assert!(
            flow.bdf1_stokes_mixed_step(&zero, &boundaries, &vec![0.0; system.n_u], 0.1,)
                .is_err()
        );
        assert_eq!(calls.get(), 0);
    }
    let expansion = |x: [f64; 2]| [x[0], 0.0];
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(e, _)| mesh.is_boundary(*e))
        .map(|(edge, _)| StokesBoundary::NavierSlip {
            edge,
            value: &expansion,
            tangential_resistance_m_s: 1.0,
        })
        .collect();
    assert!(matches!(
        system.assemble_stokes_mixed(FluxParams::default(), &zero, &boundaries),
        Err(StokesBoundaryError::IncompatibleFlux { .. })
    ));
}

/// A fixed-mesh prescribed normal flow is an explicit transport boundary.
/// Its Robin wall velocity does not replace the fluid's incoming momentum.
#[test]
fn slip_throughflow_keeps_robin_stress_separate_from_advected_momentum() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let (nu, beta, rate, crossing, lower) = (0.4, 0.8, 0.2, 0.03, 0.1);
    let tau = nu * rate;
    let upper = lower + rate + 2.0 * tau / beta;
    let exact = |x: [f64; 2]| [lower + tau / beta + rate * x[1], crossing];
    let wall = |x: [f64; 2]| [if x[1] < 0.5 { lower } else { upper }, crossing];
    let traction = |_: [f64; 2], n: [f64; 2]| [tau * n[1], tau * n[0]];
    let all: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .collect();
    let boundaries: Vec<_> = all
        .iter()
        .map(|&edge| {
            if mesh.edges[edge].normal[1].abs() > 0.5 {
                NavierStokesBoundary::NavierSlip {
                    edge,
                    value: &wall,
                    tangential_resistance_m_s: beta,
                }
            } else {
                NavierStokesBoundary::Open {
                    edge,
                    traction: &traction,
                    backflow: OpenBoundaryPolicy::NaturalTraction,
                }
            }
        })
        .collect();
    let step = system
        .bdf1_navier_stokes_mixed_step(
            FluxParams {
                nu,
                ..FluxParams::default()
            },
            &|_| [crossing * rate, 0.0],
            &boundaries,
            &affine_history(&mesh, exact),
            0.05,
            controls(),
        )
        .unwrap();
    let (error, div) = system.velocity_error(&step.solution.field.x, &exact);
    assert!(error < 3e-8 && div < 3e-8);
    let load = step
        .solution
        .reactions
        .resultant(&all, [0.0; 2], &[None])
        .unwrap();
    near(load.stress.force[0], 0.0, 6e-8);
    near(load.outward_momentum.force[0], crossing * rate, 3e-9);
    let transport = step
        .solution
        .reactions
        .transport(&all, [0.0; 2], [0.0; 2])
        .unwrap();
    near(transport.upwind_kinetic_dissipation, 0.0, 2e-9);
    near(load.slip_dissipation, 2.0 * tau * tau / beta, 3e-8);
}
