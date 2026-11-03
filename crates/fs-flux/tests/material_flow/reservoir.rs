//! Original density, viscosity and endpoint sources accompany reservoir momentum.

use super::*;

fn target(p: [f64; 2]) -> [f64; 2] {
    [0.2 * p[0] + 0.3 * p[1] + 0.4, -0.2 * p[1] + 0.1]
}

fn boundaries<'a>(
    mesh: &TriMesh,
    traction: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
) -> Vec<MaterialNavierStokesBoundary<'a>> {
    mesh.edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, geometry)| {
            if geometry.normal[0].abs() > 0.5 {
                MaterialNavierStokesBoundary::OpenReservoir {
                    edge,
                    traction,
                    velocity: &target,
                }
            } else {
                MaterialNavierStokesBoundary::Velocity {
                    edge,
                    value: &target,
                }
            }
        })
        .collect()
}

#[test]
fn actual_gas_source_scales_reservoir_momentum_stress_and_evolving_inertia_together() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let gas_spec = GasSpec::dry_air_ussa1976();
    let previous = vec![0.0; system.n_u];
    let dt = 0.025;
    let left: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, e)| mesh.is_boundary(*edge) && e.normal[0] < -0.5)
        .map(|(edge, _)| edge)
        .collect();
    for source_pressure in [101_325.0, 202_650.0] {
        let gas = GasState::try_new(&gas_spec, 293.15, source_pressure).unwrap();
        let owner =
            MaterialFluxSystem::new(&system, &gas, MaterialFluxNumerics::default()).unwrap();
        let density = owner.transport().density_kg_m3();
        let mu = owner.transport().dynamic_viscosity_pa_s();
        let force = |p: [f64; 2]| {
            let u = target(p);
            // Constant symmetric stress; grad(u)*u is exact for this affine field.
            [
                density * (u[0] / dt + 0.04 * p[0] + 0.11),
                density * (u[1] / dt + 0.04 * p[1] - 0.02),
            ]
        };
        let traction = |_: [f64; 2], n: [f64; 2]| {
            [
                (0.4 * mu - 7.0) * n[0] + 0.3 * mu * n[1],
                0.3 * mu * n[0] + (-0.4 * mu - 7.0) * n[1],
            ]
        };
        let result = owner
            .bdf1_navier_stokes_mixed_step(
                &force,
                &boundaries(&mesh, &traction),
                &previous,
                dt,
                NavierStokesControls {
                    relative_tolerance: 1e-10,
                    ..NavierStokesControls::default()
                },
            )
            .unwrap();
        let (error, divergence) = system.velocity_error(&result.solution.discrete().x, &target);
        assert!(error < 3e-9 && divergence < 3e-9, "{error}, {divergence}");
        for pressure in result.solution.pressure_pa() {
            near(*pressure, 7.0, 4e-8);
        }
        assert_eq!(
            result.solution.pressure_reference(),
            StokesPressureReference::TractionDatum
        );
        assert_eq!(result.solution.transport(), owner.transport());
        assert!(std::ptr::eq(owner.material(), &gas));
        let load = result
            .solution
            .boundary_loads(&left, [0.0; 2], &[None])
            .unwrap();
        // At x=0: u=(.4+.3y, .1-.2y). Independently integrate -rho*u_x*u
        // and its moment over y in [0,1]; the incoming normal velocity is live.
        near(
            load.outward_momentum.force_n_per_m[0],
            -0.31 * density,
            4e-8,
        );
        near(
            load.outward_momentum.force_n_per_m[1],
            0.005 * density,
            4e-8,
        );
        near(load.outward_momentum.moment_n, 73.0 * density / 400.0, 4e-8);
        near(load.stress.force_n_per_m[0], -7.0 + 0.4 * mu, 4e-8);
        near(load.stress.force_n_per_m[1], 0.3 * mu, 4e-8);
        assert_eq!(load.backflow.force_n_per_m, [0.0; 2]);
        near(
            result.solution.kinetic_energy_j_per_m(),
            density * 131.0 / 600.0,
            4e-8,
        );
        let (mass_flux, _) = result
            .solution
            .outward_mass_flux_moments(&left, [0.0; 2])
            .unwrap();
        near(mass_flux, -0.55 * density, 4e-8);
    }
    assert_eq!(previous, vec![0.0; system.n_u]);
}

#[test]
fn reservoir_rheology_uses_retained_endpoint_materials_and_preserves_source_refusals() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let previous = vec![0.0; system.n_u];
    let dt = 0.025;
    let viscosity = |rate: f64| 1.2 / (1.0 + 0.8 * rate);
    for density in [2.0, 7.0] {
        let reference = transport(density, viscosity(0.0));
        let owner = MaterialFluxSystem::new(
            &system,
            &reference,
            MaterialFluxNumerics {
                linear_tolerance: 1e-11,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        let force = |p: [f64; 2]| {
            let u = target(p);
            [
                density * (u[0] / dt + 0.04 * p[0] + 0.11),
                density * (u[1] / dt + 0.04 * p[1] - 0.02),
            ]
        };
        let traction = |_: [f64; 2], n: [f64; 2]| {
            let mu = viscosity(0.5); // sqrt(2 D:D) for the manufactured endpoint.
            [
                (0.4 * mu - 7.0) * n[0] + 0.3 * mu * n[1],
                0.3 * mu * n[0] + (-0.4 * mu - 7.0) * n[1],
            ]
        };
        let boundary = boundaries(&mesh, &traction);
        let controls = NavierStokesControls {
            relative_tolerance: 1e-10,
            relaxation: 0.8,
            ..NavierStokesControls::default()
        };
        let reads = Rc::new(Cell::new(0));
        let mut queries = 0;
        let result = owner
            .bdf1_navier_stokes_mixed_rheology_step(
                &force,
                &boundary,
                &previous,
                dt,
                controls,
                |point, retained| {
                    assert_eq!(retained, &reference);
                    queries += 1;
                    Ok(RheologyResponse {
                        transport: transport(density, viscosity(point.shear_rate_per_s)),
                        point,
                        temperature_k: 310.0,
                        reads: Rc::clone(&reads),
                    })
                },
            )
            .unwrap();
        let (error, divergence) =
            system.velocity_error(&result.step.solution.discrete().x, &target);
        assert!(error < 4e-9 && divergence < 4e-9, "{error}, {divergence}");
        assert_eq!(reads.get(), queries);
        assert!(queries > mesh.tris.len());
        for sample in result.samples() {
            assert_eq!(sample.point(), sample.material().point);
            assert_eq!(sample.transport(), &sample.material().transport);
            near(sample.point().shear_rate_per_s, 0.5, 4e-8);
            near(
                sample.transport().dynamic_viscosity_pa_s(),
                viscosity(0.5),
                4e-8,
            );
            near(sample.material().temperature_k, 310.0, 0.0);
        }
        let left: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, e)| mesh.is_boundary(*edge) && e.normal[0] < -0.5)
            .map(|(edge, _)| edge)
            .collect();
        let load = result
            .step
            .solution
            .boundary_loads(&left, [0.0; 2], &[None])
            .unwrap();
        near(
            load.outward_momentum.force_n_per_m[1],
            0.005 * density,
            4e-8,
        );
        near(load.stress.force_n_per_m[1], 0.3 * viscosity(0.5), 4e-8);
        assert_eq!(
            reads.get(),
            queries,
            "reporting must reuse accepted material owners"
        );
        let frozen = owner
            .bdf1_navier_stokes_mixed_step(&force, &boundary, &previous, dt, controls)
            .unwrap();
        let frozen_error = system
            .velocity_error(&frozen.solution.discrete().x, &target)
            .0;
        assert!(
            frozen_error > 1e-5,
            "endpoint viscosity must affect the native flow"
        );
        let refusal = MaterialError::State {
            what: "reservoir test source outside support".into(),
        };
        let failed = owner.bdf1_navier_stokes_mixed_rheology_step(
            &force,
            &boundary,
            &previous,
            dt,
            controls,
            |_, _| Err::<FluidTransportState, _>(refusal.clone()),
        );
        assert!(matches!(failed, Err(MaterialFluxError::Material(error)) if error == refusal));
        assert_eq!(previous, vec![0.0; system.n_u]);
    }
}
