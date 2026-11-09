//! G1/G3/G4: original material tangents drive pressure and slip boundaries.

use super::*;
use fs_flux::StokesRheologyControls;

fn extension_boundaries<'a>(
    mesh: &TriMesh,
    traction: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
) -> Vec<MaterialStokesBoundary<'a>> {
    mesh.edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, geometry)| {
            if geometry.normal[0] > 0.5 || geometry.normal[1] > 0.5 {
                MaterialStokesBoundary::TotalTraction {
                    edge,
                    value: traction,
                }
            } else {
                MaterialStokesBoundary::NavierSlip {
                    edge,
                    value: &no_slip,
                    tangential_resistance_pa_s_per_m: 0.0,
                }
            }
        })
        .collect()
}

#[test]
#[allow(clippy::too_many_lines)]
fn material_mixed_newton_resolves_pressure_driven_thickening_in_physical_units() {
    // Independently exact extension: u=(x,-y), gamma=2, mu=rho,
    // stress=diag(2*rho,-2*rho)-p*I. The stationary symmetry walls
    // impose no tangential velocity; the pressure difference drives the flow.
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let initial = vec![0.0; system.n_u];
    let controls = FluxNewtonControls {
        residual_tolerance: 2e-8,
        ..FluxNewtonControls::default()
    };
    for density in [2.0, 7.0] {
        let reference = transport(density, 0.2 * density);
        let owner = MaterialFluxSystem::new(
            &system,
            &reference,
            MaterialFluxNumerics {
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        for pressure_per_density in [7.0, 100_007.0] {
            let pressure = density * pressure_per_density;
            let traction = |_: [f64; 2], normal: [f64; 2]| {
                let reservoir = pressure + density * if normal[0] > 0.5 { -2.0 } else { 2.0 };
                [-reservoir * normal[0], -reservoir * normal[1]]
            };
            let boundaries = extension_boundaries(&mesh, &traction);
            if pressure_per_density == 7.0 {
                // Picard's exact scalar map a_next=5/(1+4*a^2) has
                // slope -8/5 at a=1. The same undamped budget cannot close it.
                assert!(
                    owner
                        .solve_stokes_mixed_rheology(
                            &zero_force,
                            &boundaries,
                            &initial,
                            StokesRheologyControls {
                                maximum_iterations: 12,
                                relative_tolerance: controls.residual_tolerance,
                                relaxation: 1.0,
                            },
                            |point, _| {
                                let rate = point.shear_rate_per_s;
                                Ok(transport(density, density * 0.2 * (1.0 + rate * rate)))
                            },
                        )
                        .is_err()
                );
            }
            let reads = Rc::new(Cell::new(0));
            let mut queries = 0;
            let solved = owner
                .solve_stokes_mixed_rheology_newton(
                    &zero_force,
                    &boundaries,
                    &initial,
                    controls,
                    |point, retained| {
                        assert_eq!(retained, &reference);
                        queries += 1;
                        let rate = point.shear_rate_per_s;
                        Ok(MaterialFluxViscosityTangent {
                            material: RheologyResponse {
                                transport: transport(density, density * 0.2 * (1.0 + rate * rate)),
                                point,
                                temperature_k: 315.0,
                                reads: Rc::clone(&reads),
                            },
                            derivative_pa_s2: density * 0.4 * rate,
                        })
                    },
                )
                .unwrap();
            let flow = &solved.flow;
            assert!(flow.iterations() > 1);
            if pressure_per_density == 7.0 {
                assert!(solved.backtracks > 0);
            }
            assert!(flow.nonlinear_residual() <= controls.residual_tolerance);
            assert!(flow.discrete().converged(controls.residual_tolerance));
            assert_eq!(
                flow.pressure_reference(),
                StokesPressureReference::TractionDatum
            );
            let (error, divergence) = system.velocity_error(&flow.discrete().x, &|p| [p[0], -p[1]]);
            assert!(
                error < 2e-7 && divergence < 2e-7,
                "{error:e}, {divergence:e}"
            );
            near(flow.kinetic_energy_j_per_m(), density / 3.0, density * 2e-6);
            let mut right_flux = 0.0;
            let mut top_flux = 0.0;
            for (edge, geometry) in mesh.edges.iter().enumerate() {
                if !mesh.is_boundary(edge) {
                    continue;
                }
                let flux = geometry.len * flow.velocity_dofs()[2 * edge];
                if geometry.normal[0] > 0.5 {
                    right_flux += flux;
                } else if geometry.normal[1] > 0.5 {
                    top_flux += flux;
                } else {
                    near(flux, 0.0, 1e-8);
                }
            }
            near(right_flux, 1.0, 2e-7);
            near(top_flux, -1.0, 2e-7);
            near(right_flux + top_flux, 0.0, 1e-9);
            let rates = system.shear_rates(&flow.discrete().x);
            for (cell, sample) in flow.samples().iter().enumerate() {
                assert_eq!(sample.point(), sample.material().point);
                assert_eq!(sample.point().position_m, mesh.centroids[cell]);
                assert_eq!(sample.point().shear_rate_per_s, rates[cell]);
                assert_eq!(sample.material().temperature_k, 315.0);
                assert_eq!(sample.transport().density_kg_m3(), density);
                near(
                    sample.transport().dynamic_viscosity_pa_s(),
                    density,
                    density * 2e-6,
                );
                near(
                    solved.viscosity_derivatives_pa_s2()[cell],
                    density * 0.4 * rates[cell],
                    0.0,
                );
                near(flow.pressure_pa()[cell], pressure, density * 2e-6);
                near(
                    flow.viscous_dissipation().cell_heating_w_m3()[cell],
                    4.0 * density,
                    density * 4e-6,
                );
            }
            let bulk = flow.viscous_dissipation().total_power_w_per_m();
            near(bulk, 4.0 * density, density * 4e-6);
            let pressure_work =
                -(pressure - 2.0 * density) * right_flux - (pressure + 2.0 * density) * top_flux;
            near(pressure_work, bulk, density * 5e-6);
            let right: Vec<_> = mesh
                .edges
                .iter()
                .enumerate()
                .filter(|(edge, geometry)| mesh.is_boundary(*edge) && geometry.normal[0] > 0.5)
                .map(|(edge, _)| edge)
                .collect();
            let load = flow.boundary_loads(&right, [0.0; 2], &[None]).unwrap();
            near(
                load.stress.force_n_per_m[0],
                pressure - 2.0 * density,
                density * 3e-6,
            );
            near(load.stress.force_n_per_m[1], 0.0, density * 3e-6);
            near(load.slip_dissipation_w_per_m, 0.0, 0.0);
            assert_eq!(
                reads.get(),
                queries,
                "one transport read per original owner, including reporting"
            );
        }
    }
    assert_eq!(initial, vec![0.0; system.n_u]);
}

#[test]
fn material_mixed_newton_retains_layer_sources_and_closes_bulk_and_slip_work() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let density = 4.0;
    let reference = transport(density, 0.1 * density);
    let owner = MaterialFluxSystem::new(
        &system,
        &reference,
        MaterialFluxNumerics {
            linear_tolerance: 1e-12,
            ..MaterialFluxNumerics::default()
        },
    )
    .unwrap();
    // mu=rho*K*(1+gamma^2), K=.1 below and 1.5 above y=.5.
    // gamma=3 and 1 give the SAME tau=3*rho. beta=12*rho gives slip=.25.
    let traction = |_: [f64; 2], n: [f64; 2]| {
        [
            density * (-7.0 * n[0] + 3.0 * n[1]),
            density * (3.0 * n[0] - 7.0 * n[1]),
        ]
    };
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, geometry)| {
            if geometry.normal[1] < -0.5 {
                MaterialStokesBoundary::NavierSlip {
                    edge,
                    value: &no_slip,
                    tangential_resistance_pa_s_per_m: 12.0 * density,
                }
            } else {
                MaterialStokesBoundary::TotalTraction {
                    edge,
                    value: &traction,
                }
            }
        })
        .collect();
    let reads = Rc::new(Cell::new(0));
    let mut queries = 0;
    let solved = owner
        .solve_stokes_mixed_rheology_newton(
            &zero_force,
            &boundaries,
            &vec![0.0; system.n_u],
            FluxNewtonControls::default(),
            |point, retained| {
                assert_eq!(retained, &reference);
                queries += 1;
                let lower = point.position_m[1] < 0.5;
                let coefficient = if lower { 0.1 } else { 1.5 };
                let rate = point.shear_rate_per_s;
                Ok(MaterialFluxViscosityTangent {
                    material: RheologyResponse {
                        transport: transport(density, density * coefficient * (1.0 + rate * rate)),
                        point,
                        temperature_k: if lower { 300.0 } else { 320.0 },
                        reads: Rc::clone(&reads),
                    },
                    derivative_pa_s2: density * 2.0 * coefficient * rate,
                })
            },
        )
        .unwrap();
    let flow = &solved.flow;
    let exact = |point: [f64; 2]| {
        [
            0.25 + 3.0 * point[1].min(0.5) + (point[1] - 0.5).max(0.0),
            0.0,
        ]
    };
    let (error, divergence) = system.velocity_error(&flow.discrete().x, &exact);
    assert!(
        error < 2e-7 && divergence < 2e-7,
        "{error:e}, {divergence:e}"
    );
    for (cell, sample) in flow.samples().iter().enumerate() {
        let lower = mesh.centroids[cell][1] < 0.5;
        let expected_rate = if lower { 3.0 } else { 1.0 };
        let coefficient = if lower { 0.1 } else { 1.5 };
        let rate = sample.point().shear_rate_per_s;
        assert_eq!(sample.point(), sample.material().point);
        assert_eq!(
            sample.material().temperature_k,
            if lower { 300.0 } else { 320.0 }
        );
        near(rate, expected_rate, 3e-7);
        near(
            sample.transport().dynamic_viscosity_pa_s() * rate,
            3.0 * density,
            2e-6 * density,
        );
        near(
            solved.viscosity_derivatives_pa_s2()[cell],
            density * 2.0 * coefficient * rate,
            0.0,
        );
        near(flow.pressure_pa()[cell], 7.0 * density, 2e-6 * density);
    }
    let bottom: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, geometry)| mesh.is_boundary(*edge) && geometry.normal[1] < -0.5)
        .map(|(edge, _)| edge)
        .collect();
    let top: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, geometry)| mesh.is_boundary(*edge) && geometry.normal[1] > 0.5)
        .map(|(edge, _)| edge)
        .collect();
    let wall = flow.boundary_loads(&bottom, [0.0; 2], &[None]).unwrap();
    let driving = flow.boundary_loads(&top, [0.0; 2], &[None]).unwrap();
    let bulk = flow.viscous_dissipation().total_power_w_per_m();
    near(bulk, 6.0 * density, 4e-6 * density);
    near(wall.stress.force_n_per_m[0], 3.0 * density, 2e-6 * density);
    near(
        wall.slip_dissipation_w_per_m,
        0.75 * density,
        4e-7 * density,
    );
    let input_work = -driving
        .stress
        .rigid_power_w_per_m([2.25, 0.0], 0.0)
        .unwrap();
    near(input_work, 6.75 * density, 4e-6 * density);
    near(
        input_work,
        bulk + wall.slip_dissipation_w_per_m,
        5e-6 * density,
    );
    assert_eq!(
        reads.get(),
        queries,
        "physical heat and reaction reports must not requery sources"
    );
}

#[test]
fn material_mixed_newton_preserves_thermal_carreau_source_and_analytic_derivative() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let reference = transport(2.0, 1.6);
    let owner = MaterialFluxSystem::new(
        &system,
        &reference,
        MaterialFluxNumerics {
            linear_tolerance: 1e-12,
            ..MaterialFluxNumerics::default()
        },
    )
    .unwrap();
    let law = ThermalCarreauYasuda::new(
        CarreauYasuda::new(1.6, 0.2, 3.0, 0.4, 2.0, [0.0, 1e4]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 340.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    );
    // The pressure difference is independently prescribed from the analytic
    // material at gamma=2; the solve must recover that unknown strain rate.
    let expected = law.evaluate(310.0, 2.0).unwrap();
    let exact_material = expected.bind(reference).unwrap();
    let mu = exact_material
        .transport_state()
        .unwrap()
        .dynamic_viscosity_pa_s();
    let traction = |_: [f64; 2], n: [f64; 2]| {
        let reservoir = 11.0 + if n[0] > 0.5 { -2.0 * mu } else { 2.0 * mu };
        [-reservoir * n[0], -reservoir * n[1]]
    };
    let boundaries = extension_boundaries(&mesh, &traction);
    let solved = owner
        .solve_stokes_mixed_rheology_newton(
            &zero_force,
            &boundaries,
            &vec![0.0; system.n_u],
            FluxNewtonControls::default(),
            |point, retained| {
                let state = law.evaluate(310.0, point.shear_rate_per_s)?;
                Ok(MaterialFluxViscosityTangent {
                    material: state.bind(*retained)?,
                    derivative_pa_s2: state.shear_derivative_pa_s2().unwrap(),
                })
            },
        )
        .unwrap();
    let flow = &solved.flow;
    let (error, divergence) = system.velocity_error(&flow.discrete().x, &|p| [p[0], -p[1]]);
    assert!(
        error < 2e-7 && divergence < 2e-7,
        "{error:e}, {divergence:e}"
    );
    near(
        flow.viscous_dissipation().total_power_w_per_m(),
        4.0 * mu,
        2e-7,
    );
    for (cell, sample) in flow.samples().iter().enumerate() {
        let state = sample.material().rheology();
        assert_eq!(sample.material().source(), &reference);
        assert_eq!(state.shear_rate_per_s(), sample.point().shear_rate_per_s);
        assert_eq!(state.temperature().temperature_k(), 310.0);
        assert_eq!(
            Some(solved.viscosity_derivatives_pa_s2()[cell]),
            state.shear_derivative_pa_s2()
        );
        near(flow.pressure_pa()[cell], 11.0, 2e-7);
    }
}

#[test]
fn material_mixed_newton_refuses_sources_tangents_and_density_without_publication() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let reference = transport(2.0, 0.4);
    let owner =
        MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default()).unwrap();
    let traction = |_: [f64; 2], n: [f64; 2]| {
        let reservoir = if n[0] > 0.5 { 10.0 } else { 18.0 };
        [-reservoir * n[0], -reservoir * n[1]]
    };
    let boundaries = extension_boundaries(&mesh, &traction);
    let initial = vec![0.0; system.n_u];
    let refusal = MaterialError::State {
        what: "original material refuses the pressure-driven trial strain".into(),
    };
    let mut queries = 0;
    let mut refused = false;
    let result = owner.solve_stokes_mixed_rheology_newton(
        &zero_force,
        &boundaries,
        &initial,
        FluxNewtonControls::default(),
        |point, _| {
            assert!(
                !refused,
                "source refusal must terminate even the line search"
            );
            queries += 1;
            if point.shear_rate_per_s > 0.0 {
                refused = true;
                return Err(refusal.clone());
            }
            Ok(MaterialFluxViscosityTangent {
                material: reference,
                derivative_pa_s2: 0.0,
            })
        },
    );
    assert!(matches!(result, Err(MaterialFluxError::Material(error)) if error == refusal));
    assert!(
        refused && queries > mesh.tris.len(),
        "exercise an actual Newton trial"
    );
    for (density, derivative) in [(3.0, 0.0), (2.0, f64::NAN), (2.0, 1.0)] {
        let result = owner.solve_stokes_mixed_rheology_newton(
            &zero_force,
            &boundaries,
            &initial,
            FluxNewtonControls::default(),
            |_, _| {
                Ok(MaterialFluxViscosityTangent {
                    material: transport(density, 0.4),
                    derivative_pa_s2: derivative,
                })
            },
        );
        if density != reference.density_kg_m3() {
            assert!(matches!(result, Err(MaterialFluxError::Invalid(_))));
        } else {
            assert!(matches!(
                result,
                Err(MaterialFluxError::Material(MaterialError::State { .. }))
            ));
        }
    }
    let result = owner.solve_stokes_mixed_rheology_newton(
        &zero_force,
        &boundaries,
        &initial[..system.n_u - 1],
        FluxNewtonControls::default(),
        |_, _| -> Result<MaterialFluxViscosityTangent<FluidTransportState>, MaterialError> {
            panic!("invalid initial state must precede constitutive queries")
        },
    );
    assert!(matches!(
        result,
        Err(MaterialFluxError::Newton(FluxNewtonError::Invalid(_)))
    ));
    assert_eq!(initial, vec![0.0; system.n_u]);
}
