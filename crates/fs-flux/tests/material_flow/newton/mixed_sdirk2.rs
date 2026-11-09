//! G1/G3/G4: source-owned mixed SDIRK stages, physical stress and temporal order.

use super::*;
use fs_flux::ns::newton::{FluxSdirk2Error, FluxSdirk2Stage};
use fs_flux::{MaterialNavierStokesBoundary, OpenBoundaryPolicy};

#[path = "mixed_sdirk2/adaptive.rs"]
mod adaptive;

fn moments(mesh: &TriMesh, velocity: &dyn Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    mesh.edges
        .iter()
        .flat_map(|edge| {
            let normal = |u: [f64; 2]| u[0] * edge.normal[0] + u[1] * edge.normal[1];
            let a = normal(velocity(mesh.verts[edge.verts.0]));
            let b = normal(velocity(mesh.verts[edge.verts.1]));
            [0.5 * (a + b), (b - a) / 12.0]
        })
        .collect()
}

fn stage_index(stage: FluxSdirk2Stage) -> usize {
    match stage {
        FluxSdirk2Stage::First => 0,
        FluxSdirk2Stage::Endpoint => 1,
    }
}

fn thermal_law(density: f64, layer_scale: f64) -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(
            0.2 * density * layer_scale,
            0.02 * density * layer_scale,
            0.5 * layer_scale,
            0.45,
            2.0,
            [0.0, 1e4],
        )
        .unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 350.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    )
}

fn layer_shape(y: f64) -> f64 {
    0.25 + 3.0 * y.min(0.5) + (y - 0.5).max(0.0)
}

#[test]
fn g1_mixed_sdirk2_thermal_carreau_keeps_layer_stress_slip_and_original_stages() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let time = 0.2;
    let timestep = 0.4;
    let slope = 0.4;
    let amplitude = |t: f64| 1.0 + slope * t;
    let times = [
        time + FluxSdirk2Stage::First.time_fraction() * timestep,
        time + timestep,
    ];
    let velocity = |p: [f64; 2], t: f64| [amplitude(t) * layer_shape(p[1]), 0.2];
    let previous = moments(&mesh, &|p| velocity(p, time));
    let saved = previous.clone();
    for density in [2.0, 7.0] {
        let reference = transport(density, 0.2 * density);
        let owner = MaterialFluxSystem::new(
            &system,
            &reference,
            MaterialFluxNumerics {
                penalty: 15.0,
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        let laws = [thermal_law(density, 1.0), thermal_law(density, 3.0)];
        let temperatures = times.map(|t| 300.0 + 20.0 * t);
        // Upper-layer mu(gamma) = 3*lower-layer mu(3*gamma), including
        // the original horizontal temperature shift. Both layers therefore
        // carry exactly the SAME shear stress while their gradients differ.
        let stresses: [f64; 2] = std::array::from_fn(|i| {
            let shear = 3.0 * amplitude(times[i]);
            laws[0]
                .evaluate(temperatures[i], shear)
                .unwrap()
                .dynamic_viscosity_pa_s()
                * shear
        });
        let resistance = times.map(|t| 12.0 * density * (1.0 + t));
        let pressures = times.map(|t| density * (7.0 + t));
        let tractions = [0, 1].map(|i| {
            move |_: [f64; 2], n: [f64; 2]| {
                [
                    -pressures[i] * n[0] + stresses[i] * n[1],
                    stresses[i] * n[0] - pressures[i] * n[1],
                ]
            }
        });
        let walls = [0, 1].map(|i| {
            move |p: [f64; 2]| {
                let sign = if p[1] < 0.5 { -1.0 } else { 1.0 };
                [
                    amplitude(times[i]) * layer_shape(p[1]) + sign * stresses[i] / resistance[i],
                    0.2,
                ]
            }
        });
        let reservoir = [0, 1].map(|i| move |p| velocity(p, times[i]));
        let boundaries: [Vec<_>; 2] = std::array::from_fn(|i| {
            mesh.edges
                .iter()
                .enumerate()
                .filter(|(edge, _)| mesh.is_boundary(*edge))
                .map(|(edge, geometry)| {
                    if geometry.normal[1].abs() > 0.5 {
                        MaterialNavierStokesBoundary::NavierSlip {
                            edge,
                            value: &walls[i],
                            tangential_resistance_pa_s_per_m: resistance[i],
                        }
                    } else if geometry.normal[0] < -0.5 {
                        MaterialNavierStokesBoundary::OpenReservoir {
                            edge,
                            traction: &tractions[i],
                            velocity: &reservoir[i],
                        }
                    } else {
                        MaterialNavierStokesBoundary::Open {
                            edge,
                            traction: &tractions[i],
                            backflow: OpenBoundaryPolicy::NaturalTraction,
                        }
                    }
                })
                .collect()
        });
        let force = |stage, p: [f64; 2]| {
            let i = stage_index(stage);
            let gradient = if p[1] < 0.5 { 3.0 } else { 1.0 };
            [
                density * (slope * layer_shape(p[1]) + 0.2 * amplitude(times[i]) * gradient),
                0.0,
            ]
        };
        let mut first_rates = [vec![None; mesh.tris.len()], vec![None; mesh.tris.len()]];
        let result = owner
            .sdirk2_navier_stokes_mixed_rheology_newton_step(
                &force,
                &boundaries[0],
                &boundaries[1],
                &previous,
                timestep,
                FluxNewtonControls {
                    residual_tolerance: 1e-10,
                    ..FluxNewtonControls::default()
                },
                |stage, point, retained| {
                    assert_eq!(retained, &reference);
                    let i = stage_index(stage);
                    first_rates[i][point.cell].get_or_insert(point.shear_rate_per_s);
                    let layer = usize::from(point.position_m[1] >= 0.5);
                    let state = laws[layer].evaluate(temperatures[i], point.shear_rate_per_s)?;
                    Ok(MaterialFluxViscosityTangent {
                        material: state.bind(*retained)?,
                        derivative_pa_s2: state.shear_derivative_pa_s2().unwrap(),
                    })
                },
            )
            .unwrap();
        let horizontal: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, geometry)| mesh.is_boundary(*edge) && geometry.normal[1].abs() > 0.5)
            .map(|(edge, _)| edge)
            .collect();
        let all: Vec<_> = (0..mesh.edges.len())
            .filter(|edge| mesh.is_boundary(*edge))
            .collect();
        for (i, result) in [&result.first_stage, &result.endpoint]
            .into_iter()
            .enumerate()
        {
            let flow = &result.flow.step.solution;
            let (error, divergence) =
                system.velocity_error(&flow.discrete().x, &|p| velocity(p, times[i]));
            assert!(
                error < 4e-7 && divergence < 4e-7,
                "{error:e}, {divergence:e}"
            );
            assert_eq!(
                flow.pressure_reference(),
                StokesPressureReference::TractionDatum
            );
            near(
                flow.viscous_dissipation().total_power_w_per_m(),
                2.0 * amplitude(times[i]) * stresses[i],
                2e-5 * density,
            );
            let loads = flow.boundary_loads(&horizontal, [0.0; 2], &[None]).unwrap();
            near(
                loads.slip_dissipation_w_per_m,
                2.0 * stresses[i].powi(2) / resistance[i],
                2e-5 * density,
            );
            let exchange = flow.boundary_transport(&all, [0.0; 2], [0.0; 2]).unwrap();
            near(
                exchange.outward_kinetic_energy_w_per_m,
                0.5 * density * amplitude(times[i]).powi(2),
                2e-5 * density,
            );
            near(
                flow.kinetic_energy_j_per_m(),
                0.5 * density * (125.0 / 48.0 * amplitude(times[i]).powi(2) + 0.04),
                2e-5 * density,
            );
            for (cell, sample) in result.flow.samples().iter().enumerate() {
                let gradient = if mesh.centroids[cell][1] < 0.5 {
                    3.0
                } else {
                    1.0
                };
                let initial_time = if i == 0 { time } else { times[0] };
                near(
                    first_rates[i][cell].unwrap(),
                    gradient * amplitude(initial_time),
                    2e-6,
                );
                let state = sample.material().rheology();
                assert_eq!(sample.material().source(), &reference);
                assert_eq!(state.temperature().temperature_k(), temperatures[i]);
                assert_eq!(state.shear_rate_per_s(), sample.point().shear_rate_per_s);
                assert_eq!(
                    state.shear_derivative_pa_s2(),
                    Some(result.viscosity_derivatives_pa_s2()[cell])
                );
                near(
                    sample.point().shear_rate_per_s,
                    gradient * amplitude(times[i]),
                    2e-6,
                );
                near(
                    sample.transport().dynamic_viscosity_pa_s() * gradient * amplitude(times[i]),
                    stresses[i],
                    2e-5 * density,
                );
                near(flow.pressure_pa()[cell], pressures[i], 2e-5 * density);
            }
        }
    }
    assert_eq!(previous, saved);
}

#[test]
fn g1_mixed_sdirk2_has_second_order_time_error_without_prescribed_velocity() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let density = 3.0;
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
    let law = thermal_law(density, 1.0);
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|edge| mesh.is_boundary(*edge))
        .map(|edge| MaterialNavierStokesBoundary::Open {
            edge,
            traction: &traction,
            backflow: OpenBoundaryPolicy::NaturalTraction,
        })
        .collect();
    let final_time = 0.4;
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let dt = final_time / steps as f64;
        let mut previous = moments(&mesh, &|_| [1.0, 0.0]);
        for step in 0..steps {
            let time = step as f64 * dt;
            let result = owner
                .sdirk2_navier_stokes_mixed_rheology_newton_step(
                    &|stage, _| {
                        [
                            3.0 * density * (time + stage.time_fraction() * dt).powi(2),
                            0.0,
                        ]
                    },
                    &boundaries,
                    &boundaries,
                    &previous,
                    dt,
                    FluxNewtonControls {
                        residual_tolerance: 1e-10,
                        ..FluxNewtonControls::default()
                    },
                    |stage, point, retained| {
                        let temperature = 300.0 + 20.0 * (time + stage.time_fraction() * dt);
                        let state = law.evaluate(temperature, point.shear_rate_per_s)?;
                        Ok(MaterialFluxViscosityTangent {
                            material: state.bind(*retained)?,
                            derivative_pa_s2: state.shear_derivative_pa_s2().unwrap(),
                        })
                    },
                )
                .unwrap();
            previous = result.endpoint.flow.step.solution.discrete().x[..system.n_u].to_vec();
        }
        // The exact solution is uniform translation u=1+t³. Stress and
        // convection have zero divergence; no prescribed velocity can mask
        // an erroneous inertial history or stage timestamp.
        let (error, divergence) =
            system.velocity_error(&previous, &|_| [1.0 + final_time.powi(3), 0.0]);
        assert!(divergence < 1e-8);
        let analytic_error = (1.5 * core::f64::consts::FRAC_1_SQRT_2 - 1.0) * final_time * dt * dt;
        near(error, analytic_error, 2e-8);
        errors.push(error);
    }
    for pair in errors.windows(2) {
        assert!((pair[0] / pair[1] - 4.0).abs() < 0.02, "{errors:?}");
    }
}

#[test]
fn g4_mixed_sdirk2_preserves_stage_refusals_history_and_clean_retry() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let reference = transport(2.0, 0.4);
    let owner =
        MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default()).unwrap();
    let previous = vec![0.0; system.n_u];
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|edge| mesh.is_boundary(*edge))
        .map(|edge| MaterialNavierStokesBoundary::Open {
            edge,
            traction: &traction,
            backflow: OpenBoundaryPolicy::NaturalTraction,
        })
        .collect();
    let refusal = MaterialError::State {
        what: "endpoint temperature outside original source".into(),
    };
    let mut queries = [0; 2];
    let result = owner.sdirk2_navier_stokes_mixed_rheology_newton_step(
        &|_, _| [0.2, 0.0],
        &boundaries,
        &boundaries,
        &previous,
        0.1,
        FluxNewtonControls::default(),
        |stage, _, _| {
            queries[stage_index(stage)] += 1;
            if stage == FluxSdirk2Stage::Endpoint {
                return Err(refusal.clone());
            }
            Ok(MaterialFluxViscosityTangent {
                material: reference,
                derivative_pa_s2: 0.0,
            })
        },
    );
    assert!(queries[0] >= mesh.tris.len());
    assert_eq!(queries[1], 1, "a source refusal is terminal");
    assert!(
        matches!(result, Err(FluxSdirk2Error { stage: Some(FluxSdirk2Stage::Endpoint), error: MaterialFluxError::Material(error) }) if error == refusal)
    );
    assert!(previous.iter().all(|value| *value == 0.0));
    let invalid = owner.sdirk2_navier_stokes_mixed_rheology_newton_step(
        &|_, _| panic!("invalid timestep queried a force"),
        &boundaries,
        &boundaries,
        &previous,
        0.0,
        FluxNewtonControls::default(),
        |_, _, _| -> Result<MaterialFluxViscosityTangent<FluidTransportState>, MaterialError> {
            panic!("invalid timestep queried a source")
        },
    );
    assert!(matches!(
        invalid,
        Err(FluxSdirk2Error {
            stage: None,
            error: MaterialFluxError::Newton(FluxNewtonError::Invalid(_))
        })
    ));
    let bad_boundary = owner.sdirk2_navier_stokes_mixed_rheology_newton_step(
        &|_, _| [0.2, 0.0],
        &boundaries,
        &[],
        &previous,
        0.1,
        FluxNewtonControls::default(),
        |_, _, _| {
            Ok(MaterialFluxViscosityTangent {
                material: reference,
                derivative_pa_s2: 0.0,
            })
        },
    );
    assert!(matches!(
        bad_boundary,
        Err(FluxSdirk2Error {
            stage: Some(FluxSdirk2Stage::Endpoint),
            error: MaterialFluxError::Newton(FluxNewtonError::Boundary(_))
        })
    ));
    let retry = owner
        .sdirk2_navier_stokes_mixed_rheology_newton_step(
            &|_, _| [0.2, 0.0],
            &boundaries,
            &boundaries,
            &previous,
            0.1,
            FluxNewtonControls::default(),
            |_, _, _| {
                Ok(MaterialFluxViscosityTangent {
                    material: reference,
                    derivative_pa_s2: 0.0,
                })
            },
        )
        .unwrap();
    let (error, divergence) = system
        .velocity_error(&retry.endpoint.flow.step.solution.discrete().x, &|_| {
            [0.01, 0.0]
        });
    assert!(error < 1e-8 && divergence < 1e-8);
    assert!(previous.iter().all(|value| *value == 0.0));
}
