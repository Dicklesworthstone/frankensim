//! G1/G3/G4: original material tangents, mixed transport and physical inertia.

use super::*;
use fs_flux::{MaterialNavierStokesBoundary, NavierStokesControls, OpenBoundaryPolicy};

fn moments(mesh: &TriMesh, velocity: &dyn Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    mesh.edges
        .iter()
        .flat_map(|edge| {
            let a = velocity(mesh.verts[edge.verts.0]);
            let b = velocity(mesh.verts[edge.verts.1]);
            let normal = |u: [f64; 2]| u[0] * edge.normal[0] + u[1] * edge.normal[1];
            // Exact mean and centered first moment of each affine normal trace.
            [
                0.5 * (normal(a) + normal(b)),
                (normal(b) - normal(a)) / 12.0,
            ]
        })
        .collect()
}

fn layered_velocity(p: [f64; 2], amplitude: f64) -> [f64; 2] {
    [
        amplitude * (0.25 + 3.0 * p[1].min(0.5) + (p[1] - 0.5).max(0.0)),
        0.2,
    ]
}

#[test]
fn g1_material_mixed_newton_bdf_keeps_layer_stress_slip_transport_and_history() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let previous = moments(&mesh, &|p| layered_velocity(p, 0.8_f64.powi(2)));
    let older = moments(&mesh, &|p| layered_velocity(p, 0.5_f64.powi(2)));
    let saved = (previous.clone(), older.clone());
    let target = |p| layered_velocity(p, 1.0);
    let walls = |p: [f64; 2]| [if p[1] < 0.5 { 0.0 } else { 2.5 }, 0.2];
    let controls = FluxNewtonControls {
        residual_tolerance: 1e-10,
        ..FluxNewtonControls::default()
    };
    for density in [2.0, 7.0] {
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
                if geometry.normal[1].abs() > 0.5 {
                    // Prescribed-normal transport, not an impermeable wall.
                    MaterialNavierStokesBoundary::NavierSlip {
                        edge,
                        value: &walls,
                        tangential_resistance_pa_s_per_m: 12.0 * density,
                    }
                } else if geometry.normal[0] < -0.5 {
                    MaterialNavierStokesBoundary::OpenReservoir {
                        edge,
                        traction: &traction,
                        velocity: &target,
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
        for second_order in [false, true] {
            // Unequal BDF2 is exact for a(t)=t² at t=1, .8 and .5.
            // BDF1 uses its own actual difference, with no fictitious history.
            let acceleration = if second_order {
                2.0
            } else {
                (1.0 - 0.8_f64.powi(2)) / 0.2
            };
            let force = |p: [f64; 2]| {
                let slope = if p[1] < 0.5 { 3.0 } else { 1.0 };
                [density * (acceleration * target(p)[0] + 0.2 * slope), 0.0]
            };
            let reads = Rc::new(Cell::new(0));
            let mut queries = 0;
            let mut first_rates = vec![None; mesh.tris.len()];
            let result = owner
                .bdf2_navier_stokes_mixed_rheology_newton_variable_step(
                    &force,
                    &boundaries,
                    &previous,
                    second_order.then_some((older.as_slice(), 0.3)),
                    0.2,
                    controls,
                    |point, retained| {
                        assert_eq!(retained, &reference);
                        first_rates[point.cell].get_or_insert(point.shear_rate_per_s);
                        queries += 1;
                        let lower = point.position_m[1] < 0.5;
                        let coefficient = if lower { 0.1 } else { 1.5 };
                        let rate = point.shear_rate_per_s;
                        Ok(MaterialFluxViscosityTangent {
                            material: RheologyResponse {
                                transport: transport(
                                    density,
                                    density * coefficient * (1.0 + rate * rate),
                                ),
                                point,
                                temperature_k: if lower { 300.0 } else { 320.0 },
                                reads: Rc::clone(&reads),
                            },
                            derivative_pa_s2: 2.0 * density * coefficient * rate,
                        })
                    },
                )
                .unwrap();
            let flow = &result.flow.step.solution;
            let (error, divergence) = system.velocity_error(&flow.discrete().x, &target);
            assert!(
                error < 3e-7 && divergence < 3e-7,
                "{error:e}, {divergence:e}"
            );
            assert!(result.flow.step.iterations > 1);
            assert!(result.flow.step.nonlinear_residual <= controls.residual_tolerance);
            assert_eq!(
                flow.pressure_reference(),
                StokesPressureReference::TractionDatum
            );
            // Reservoir inflow keeps literal stress; its policy is NaturalTraction.
            assert!(
                result
                    .flow
                    .step
                    .open_boundary_policies
                    .iter()
                    .all(|(_, p)| *p == OpenBoundaryPolicy::NaturalTraction)
            );
            let rates = system.shear_rates(&flow.discrete().x);
            for (cell, sample) in result.flow.samples().iter().enumerate() {
                let lower = mesh.centroids[cell][1] < 0.5;
                let coefficient = if lower { 0.1 } else { 1.5 };
                let shear = if lower { 3.0 } else { 1.0 };
                near(first_rates[cell].unwrap(), 0.8_f64.powi(2) * shear, 2e-12);
                assert_eq!(sample.point(), sample.material().point);
                assert_eq!(sample.point().position_m, mesh.centroids[cell]);
                assert_eq!(sample.point().shear_rate_per_s, rates[cell]);
                assert_eq!(
                    sample.material().temperature_k,
                    if lower { 300.0 } else { 320.0 }
                );
                near(sample.point().shear_rate_per_s, shear, 8e-7);
                near(
                    sample.transport().dynamic_viscosity_pa_s() * shear,
                    3.0 * density,
                    1e-5 * density,
                );
                near(
                    result.viscosity_derivatives_pa_s2()[cell],
                    2.0 * density * coefficient * rates[cell],
                    0.0,
                );
                near(flow.pressure_pa()[cell], 7.0 * density, 1e-5 * density);
            }
            let all: Vec<_> = (0..mesh.edges.len())
                .filter(|e| mesh.is_boundary(*e))
                .collect();
            let bottom: Vec<_> = all
                .iter()
                .copied()
                .filter(|e| mesh.edges[*e].normal[1] < -0.5)
                .collect();
            let top: Vec<_> = all
                .iter()
                .copied()
                .filter(|e| mesh.edges[*e].normal[1] > 0.5)
                .collect();
            let low = flow.boundary_loads(&bottom, [0.0; 2], &[None]).unwrap();
            let high = flow.boundary_loads(&top, [0.0; 2], &[None]).unwrap();
            let exchange = flow.boundary_transport(&all, [0.0; 2], [0.0; 2]).unwrap();
            let bulk = flow.viscous_dissipation().total_power_w_per_m();
            let slip = low.slip_dissipation_w_per_m + high.slip_dissipation_w_per_m;
            near(bulk, 6.0 * density, 1e-5 * density);
            near(slip, 1.5 * density, 1e-5 * density);
            near(low.stress.force_n_per_m[0], 3.0 * density, 1e-5 * density);
            near(high.stress.force_n_per_m[0], -3.0 * density, 1e-5 * density);
            near(exchange.outward_mass_flux_kg_per_m_s, 0.0, 1e-6 * density);
            near(
                exchange.outward_kinetic_energy_w_per_m,
                0.5 * density,
                1e-5 * density,
            );
            assert!(exchange.upwind_kinetic_dissipation_w_per_m < 1e-10 * density);
            // Integral h² dy = 125/48, independently of the discrete solver.
            let shape_energy = 125.0 / 48.0;
            near(
                flow.kinetic_energy_j_per_m(),
                0.5 * density * (shape_energy + 0.04),
                1e-5 * density,
            );
            if !second_order {
                let before = 0.5 * density * (0.8_f64.powi(4) * shape_energy + 0.04);
                let bdf_damping = 0.5 * density * (1.0 - 0.8_f64.powi(2)).powi(2) * shape_energy;
                let body_work = density * (acceleration * shape_energy + 0.5);
                let wall_work = 7.5 * density;
                let balance = flow.kinetic_energy_j_per_m() - before
                    + bdf_damping
                    + 0.2
                        * (bulk + slip + exchange.outward_kinetic_energy_w_per_m
                            - body_work
                            - wall_work);
                near(balance, 0.0, 2e-5 * density);
            }
            assert_eq!(
                reads.get(),
                queries,
                "reporting must retain each original source response"
            );
        }
    }
    assert_eq!((previous, older), saved);
}

#[test]
fn g1_original_thermal_carreau_newton_accelerates_an_all_open_component() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let previous = moments(&mesh, &|p| [0.2 * p[0], -0.2 * p[1]]);
    let saved = previous.clone();
    let target = |p: [f64; 2]| [0.4 * p[0], -0.4 * p[1]];
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
    let expected = law.evaluate(310.0, 0.8).unwrap();
    let mu = expected.dynamic_viscosity_pa_s();
    for density in [2.0, 7.0] {
        let reference = transport(density, 1.6);
        let owner = MaterialFluxSystem::new(
            &system,
            &reference,
            MaterialFluxNumerics {
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        let force = |p: [f64; 2]| {
            [
                density * (2.0 + 0.16) * p[0],
                density * (-2.0 + 0.16) * p[1],
            ]
        };
        let traction =
            |_: [f64; 2], n: [f64; 2]| [(0.8 * mu - 11.0) * n[0], (-0.8 * mu - 11.0) * n[1]];
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, _)| MaterialNavierStokesBoundary::OpenReservoir {
                edge,
                traction: &traction,
                velocity: &target,
            })
            .collect();
        let result = owner
            .bdf1_navier_stokes_mixed_rheology_newton_step(
                &force,
                &boundaries,
                &previous,
                0.1,
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
        let flow = &result.flow.step.solution;
        let (error, divergence) = system.velocity_error(&flow.discrete().x, &target);
        assert!(
            error < 2e-7 && divergence < 2e-7,
            "{error:e}, {divergence:e}"
        );
        near(
            flow.viscous_dissipation().total_power_w_per_m(),
            0.64 * mu,
            2e-7,
        );
        near(
            flow.kinetic_energy_j_per_m(),
            density * 0.16 / 3.0,
            2e-7 * density,
        );
        for (cell, sample) in result.flow.samples().iter().enumerate() {
            let state = sample.material().rheology();
            assert_eq!(sample.material().source(), &reference);
            assert_eq!(state.temperature().temperature_k(), 310.0);
            assert_eq!(state.shear_rate_per_s(), sample.point().shear_rate_per_s);
            assert_eq!(
                state.shear_derivative_pa_s2(),
                Some(result.viscosity_derivatives_pa_s2()[cell])
            );
            near(sample.point().shear_rate_per_s, 0.8, 2e-7);
            near(flow.pressure_pa()[cell], 11.0, 2e-7);
        }
        assert!(
            result
                .flow
                .step
                .open_boundary_policies
                .iter()
                .all(|(_, p)| *p == OpenBoundaryPolicy::NaturalTraction)
        );
    }
    assert_eq!(previous, saved);
}

#[test]
fn g3_material_steady_mixed_newton_preserves_live_convection_and_backflow_policies() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let target = |p: [f64; 2]| [0.4 + 0.3 * p[1], 0.1];
    let initial = vec![0.0; system.n_u];
    let density = 3.0;
    let mu = density * (0.2 + 0.3 * 0.3_f64.powi(2));
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
    for policy in [
        OpenBoundaryPolicy::NaturalTraction,
        OpenBoundaryPolicy::DissipativeBackflow,
    ] {
        let traction = |p: [f64; 2], n: [f64; 2]| {
            let u = target(p);
            let inward = (u[0] * n[0] + u[1] * n[1]).min(0.0);
            let undo_backflow = if policy == OpenBoundaryPolicy::DissipativeBackflow {
                density * inward
            } else {
                0.0
            };
            [
                -7.0 * density * n[0] + 0.3 * mu * n[1] - undo_backflow * u[0],
                0.3 * mu * n[0] - 7.0 * density * n[1] - undo_backflow * u[1],
            ]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, geometry)| {
                if geometry.normal[1].abs() > 0.5 {
                    MaterialNavierStokesBoundary::Velocity {
                        edge,
                        value: &target,
                    }
                } else {
                    MaterialNavierStokesBoundary::Open {
                        edge,
                        traction: &traction,
                        backflow: policy,
                    }
                }
            })
            .collect();
        let force = |_: [f64; 2]| [0.03 * density, 0.0];
        let result = owner
            .solve_navier_stokes_mixed_rheology_newton(
                &force,
                &boundaries,
                &initial,
                FluxNewtonControls::default(),
                |point, _| {
                    let rate = point.shear_rate_per_s;
                    Ok(MaterialFluxViscosityTangent {
                        material: transport(density, density * (0.2 + 0.3 * rate * rate)),
                        derivative_pa_s2: density * 0.6 * rate,
                    })
                },
            )
            .unwrap();
        let picard = owner
            .solve_navier_stokes_mixed_rheology(
                &force,
                &boundaries,
                &initial,
                NavierStokesControls {
                    maximum_iterations: 100,
                    relative_tolerance: 1e-10,
                    relaxation: 1.0,
                },
                |point, _| {
                    Ok(transport(
                        density,
                        density * (0.2 + 0.3 * point.shear_rate_per_s.powi(2)),
                    ))
                },
            )
            .unwrap();
        let flow = &result.flow.step.solution;
        let (error, divergence) = system.velocity_error(&flow.discrete().x, &target);
        assert!(
            error < 2e-7 && divergence < 2e-7,
            "{error:e}, {divergence:e}"
        );
        for (actual, reference) in flow
            .discrete()
            .x
            .iter()
            .zip(&picard.step.solution.discrete().x)
        {
            near(*actual, *reference, 2e-7);
        }
        near(
            flow.kinetic_energy_j_per_m(),
            0.16 * density,
            2e-7 * density,
        );
        near(
            flow.viscous_dissipation().total_power_w_per_m(),
            0.09 * mu,
            2e-7,
        );
        assert!(
            result
                .flow
                .step
                .open_boundary_policies
                .iter()
                .all(|(_, p)| *p == policy)
        );
        for (cell, sample) in result.flow.samples().iter().enumerate() {
            near(flow.pressure_pa()[cell], 7.0 * density, 2e-7 * density);
            near(
                result.viscosity_derivatives_pa_s2()[cell],
                density * 0.6 * sample.point().shear_rate_per_s,
                0.0,
            );
        }
    }
    assert_eq!(initial, vec![0.0; system.n_u]);
}

#[test]
fn g4_material_mixed_newton_preserves_original_refusals_and_all_histories() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let previous = vec![0.0; system.n_u];
    let older = previous.clone();
    let density = 2.0;
    let reference = transport(density, 0.2 * density);
    let owner =
        MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default()).unwrap();
    let target = |p: [f64; 2]| [0.4 + 0.3 * p[1], 0.1];
    let calls = Cell::new(0);
    let mu = density * (0.2 + 0.8 * 0.3_f64.powi(2));
    let traction = |_: [f64; 2], n: [f64; 2]| {
        calls.set(calls.get() + 1);
        [
            -7.0 * density * n[0] + 0.3 * mu * n[1],
            0.3 * mu * n[0] - 7.0 * density * n[1],
        ]
    };
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, geometry)| {
            if geometry.normal[1].abs() > 0.5 {
                MaterialNavierStokesBoundary::Velocity {
                    edge,
                    value: &target,
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
    let force = |p: [f64; 2]| {
        calls.set(calls.get() + 1);
        [density * (target(p)[0] / 0.1 + 0.03), density]
    };
    for (history, interval) in [
        (Some((&older[..1], 0.1)), 0.1),
        (Some((older.as_slice(), 0.1)), 0.21),
        (Some((older.as_slice(), f64::NAN)), 0.1),
        (None, 0.0),
    ] {
        let result = owner.bdf2_navier_stokes_mixed_rheology_newton_variable_step(
            &force,
            &boundaries,
            &previous,
            history,
            interval,
            FluxNewtonControls::default(),
            |_, _| {
                calls.set(calls.get() + 1);
                Ok(MaterialFluxViscosityTangent {
                    material: reference,
                    derivative_pa_s2: 0.0,
                })
            },
        );
        assert!(matches!(result, Err(MaterialFluxError::Newton(_))));
        assert_eq!(
            calls.get(),
            0,
            "invalid history must precede original callbacks"
        );
    }
    let refusal = MaterialError::State {
        what: "original mixed-flow rate support".into(),
    };
    let mut queries = 0;
    let result = owner.bdf1_navier_stokes_mixed_rheology_newton_step(
        &force,
        &boundaries,
        &previous,
        0.1,
        FluxNewtonControls::default(),
        |point, _| {
            queries += 1;
            if point.shear_rate_per_s > 1e-12 {
                return Err(refusal.clone());
            }
            Ok(MaterialFluxViscosityTangent {
                material: reference,
                derivative_pa_s2: 0.0,
            })
        },
    );
    assert!(
        queries > mesh.tris.len(),
        "exercise a real Newton trial after initial source admission"
    );
    assert!(matches!(result, Err(MaterialFluxError::Material(error)) if error == refusal));
    for invalid_derivative in [false, true] {
        let result = owner.bdf1_navier_stokes_mixed_rheology_newton_step(
            &force,
            &boundaries,
            &previous,
            0.1,
            FluxNewtonControls::default(),
            |_, _| {
                Ok(MaterialFluxViscosityTangent {
                    material: if invalid_derivative {
                        reference
                    } else {
                        transport(3.0, 0.2)
                    },
                    derivative_pa_s2: if invalid_derivative { f64::NAN } else { 0.0 },
                })
            },
        );
        assert!(matches!(
            result,
            Err(MaterialFluxError::Material(_) | MaterialFluxError::Invalid(_))
        ));
    }
    let source = |point: MaterialFluxRheologyPoint, _: &FluidTransportState| {
        let rate = point.shear_rate_per_s;
        Ok(MaterialFluxViscosityTangent {
            material: transport(density, density * (0.2 + 0.8 * rate * rate)),
            derivative_pa_s2: 1.6 * density * rate,
        })
    };
    let retry = owner
        .bdf1_navier_stokes_mixed_rheology_newton_step(
            &force,
            &boundaries,
            &previous,
            0.1,
            FluxNewtonControls::default(),
            source,
        )
        .unwrap();
    let startup = owner
        .bdf2_navier_stokes_mixed_rheology_newton_variable_step(
            &force,
            &boundaries,
            &previous,
            None,
            0.1,
            FluxNewtonControls::default(),
            source,
        )
        .unwrap();
    assert_eq!(
        retry.flow.step.solution.discrete().x,
        startup.flow.step.solution.discrete().x
    );
    assert_eq!(
        retry.viscosity_derivatives_pa_s2(),
        startup.viscosity_derivatives_pa_s2()
    );
    assert_eq!(previous, vec![0.0; system.n_u]);
    assert_eq!(older, previous);
}
