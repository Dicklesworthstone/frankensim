//! G1/G3/G4: adaptive physical-time ownership, mass norms and atomic refusals.

use super::*;
use fs_exec::CancelGate;
use fs_flux::material::adaptive::*;

#[path = "adaptive/heat_control.rs"]
mod heat_control;

fn controls(atol: f64) -> MaterialSdirk2AdaptiveControls {
    MaterialSdirk2AdaptiveControls {
        initial_step_s: 0.4,
        minimum_step_s: 1e-8,
        maximum_step_s: 0.4,
        velocity_absolute_tolerance_m_s: atol,
        velocity_relative_tolerance: 0.0,
        maximum_attempts: 500,
        newton: FluxNewtonControls {
            residual_tolerance: 1e-11,
            ..FluxNewtonControls::default()
        },
        ..MaterialSdirk2AdaptiveControls::default()
    }
}

#[test]
fn g1_adaptive_original_thermal_carreau_uses_physical_rms_and_absolute_stage_times() {
    let gate = CancelGate::new_clock_free();
    let mut fine_reports = Vec::new();
    for (width, height, nx, ny) in [(1.0, 1.0, 1, 1), (2.0, 3.0, 2, 3)] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(width, height, nx, ny));
        let system = FluxSystem::new(&mesh);
        let reference = transport(3.0, 0.6);
        let owner = MaterialFluxSystem::new(
            &system,
            &reference,
            MaterialFluxNumerics {
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        let law = thermal_law(3.0, 1.0);
        let start = 2.0;
        let end = 2.4;
        let previous = moments(&mesh, &|_| [1.0, 0.0]);
        let saved = previous.clone();
        let traction = |_: f64, _: [f64; 2], _: [f64; 2]| [0.0; 2];
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|edge| mesh.is_boundary(*edge))
            .map(|edge| MaterialNavierStokesTimeBoundary::Open {
                edge,
                traction: &traction,
                backflow: OpenBoundaryPolicy::NaturalTraction,
            })
            .collect();
        let mut errors = Vec::new();
        for tolerance in [1e-4, 1e-6] {
            let result = owner
                .advance_navier_stokes_mixed_rheology_sdirk2(
                    start,
                    end,
                    &previous,
                    &|t, _| [9.0 * (t - start).powi(2), 0.0],
                    &boundaries,
                    &controls(tolerance),
                    |t, point, retained| {
                        assert!(t > start && t <= end);
                        assert_eq!(retained, &reference);
                        let state =
                            law.evaluate(300.0 + 10.0 * (t - start), point.shear_rate_per_s)?;
                        Ok(MaterialFluxViscosityTangent {
                            material: state.bind(*retained)?,
                            derivative_pa_s2: state.shear_derivative_pa_s2().unwrap(),
                        })
                    },
                    &gate,
                )
                .unwrap();
            let report = result.report;
            assert_eq!(report.start_time_s, start);
            assert_eq!(report.end_time_s, end);
            assert!(report.accepted_comparisons > 1 && report.rejected_comparisons > 0);
            assert_eq!(report.numerical_rejections, 0);
            assert_eq!(report.stage_solves, 6 * report.attempts);
            assert!(report.maximum_accepted_error <= 1.0);
            let (l2, divergence) = system.velocity_error(
                &result.second_half.endpoint.flow.step.solution.discrete().x,
                &|_| [1.0 + (end - start).powi(3), 0.0],
            );
            let rms = l2 / (width * height).sqrt();
            errors.push(rms);
            assert!(divergence < 1e-8);
            let mid = report.last_step_start_time_s + 0.5 * (end - report.last_step_start_time_s);
            let first_stage_time = |a, b| a + FluxSdirk2Stage::First.time_fraction() * (b - a);
            for (time, stage) in [
                (
                    first_stage_time(report.last_step_start_time_s, mid),
                    &result.first_half.first_stage,
                ),
                (mid, &result.first_half.endpoint),
                (first_stage_time(mid, end), &result.second_half.first_stage),
                (end, &result.second_half.endpoint),
            ] {
                for sample in stage.flow.samples() {
                    assert_eq!(sample.material().source(), &reference);
                    near(
                        sample.material().rheology().temperature().temperature_k(),
                        300.0 + 10.0 * (time - start),
                        1e-12,
                    );
                }
            }
            if tolerance == 1e-6 {
                fine_reports.push((report, rms));
            }
        }
        // Exact cubic translation is a spatially represented solution on both
        // meshes; tighter LOCAL tolerance reduces its accumulated time error.
        assert!(
            errors[1] < 1e-5 && errors[1] < 0.1 * errors[0],
            "{errors:?}"
        );
        assert_eq!(previous, saved);
    }
    // Changing domain area and edge-moment count must not rescale the physical
    // RMS tolerance. Raw coefficient norms would produce a different controller.
    assert_eq!(
        fine_reports[0].0.accepted_comparisons,
        fine_reports[1].0.accepted_comparisons
    );
    assert_eq!(
        fine_reports[0].0.rejected_comparisons,
        fine_reports[1].0.rejected_comparisons
    );
    near(fine_reports[0].1, fine_reports[1].1, 2e-8);
}

#[test]
fn g1_adaptive_mixed_carreau_retains_time_dependent_slip_traction_reservoir_and_sources() {
    let gate = CancelGate::new_clock_free();
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let density = 2.0;
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
    let amplitude = |time: f64| 1.0 + 0.4 * (time - 1.0);
    let temperature = |time: f64| 300.0 + 20.0 * (time - 1.0);
    let velocity = |t: f64, p: [f64; 2]| [amplitude(t) * layer_shape(p[1]), 0.2];
    let pressure = |time: f64| density * (7.0 + time);
    let stress = |time: f64| {
        let rate = 3.0 * amplitude(time);
        laws[0]
            .evaluate(temperature(time), rate)
            .unwrap()
            .dynamic_viscosity_pa_s()
            * rate
    };
    let resistance = |time| 12.0 * density * (1.0 + time);
    let traction = |time, _: [f64; 2], n: [f64; 2]| {
        [
            -pressure(time) * n[0] + stress(time) * n[1],
            stress(time) * n[0] - pressure(time) * n[1],
        ]
    };
    let wall = |time, p: [f64; 2]| {
        let sign = if p[1] < 0.5 { -1.0 } else { 1.0 };
        [
            velocity(time, p)[0] + sign * stress(time) / resistance(time),
            0.2,
        ]
    };
    let boundaries: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, geometry)| {
            if geometry.normal[1].abs() > 0.5 {
                MaterialNavierStokesTimeBoundary::NavierSlip {
                    edge,
                    value: &wall,
                    tangential_resistance_pa_s_per_m: &resistance,
                }
            } else if geometry.normal[0] < -0.5 {
                MaterialNavierStokesTimeBoundary::OpenReservoir {
                    edge,
                    traction: &traction,
                    velocity: &velocity,
                }
            } else {
                MaterialNavierStokesTimeBoundary::Open {
                    edge,
                    traction: &traction,
                    backflow: OpenBoundaryPolicy::NaturalTraction,
                }
            }
        })
        .collect();
    let previous = moments(&mesh, &|p| velocity(1.0, p));
    let saved = previous.clone();
    let mut config = controls(1e-5);
    config.initial_step_s = 0.2;
    config.maximum_step_s = 0.2;
    let result = owner
        .advance_navier_stokes_mixed_rheology_sdirk2(
            1.0,
            1.4,
            &previous,
            &|time, p| {
                let gradient = if p[1] < 0.5 { 3.0 } else { 1.0 };
                [
                    density * (0.4 * layer_shape(p[1]) + 0.2 * amplitude(time) * gradient),
                    0.0,
                ]
            },
            &boundaries,
            &config,
            |time, point, retained| {
                let layer = usize::from(point.position_m[1] >= 0.5);
                let state = laws[layer].evaluate(temperature(time), point.shear_rate_per_s)?;
                Ok(MaterialFluxViscosityTangent {
                    material: state.bind(*retained)?,
                    derivative_pa_s2: state.shear_derivative_pa_s2().unwrap(),
                })
            },
            &gate,
        )
        .unwrap();
    let endpoint = &result.second_half.endpoint;
    let flow = &endpoint.flow.step.solution;
    let (error, divergence) = system.velocity_error(&flow.discrete().x, &|p| velocity(1.4, p));
    assert!(error < 2e-6 && divergence < 2e-6);
    near(
        flow.viscous_dissipation().total_power_w_per_m(),
        2.0 * amplitude(1.4) * stress(1.4),
        4e-5,
    );
    let horizontal: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, e)| mesh.is_boundary(*edge) && e.normal[1].abs() > 0.5)
        .map(|(edge, _)| edge)
        .collect();
    near(
        flow.boundary_loads(&horizontal, [0.0; 2], &[None])
            .unwrap()
            .slip_dissipation_w_per_m,
        2.0 * stress(1.4).powi(2) / resistance(1.4),
        4e-5,
    );
    for (cell, sample) in endpoint.flow.samples().iter().enumerate() {
        assert_eq!(sample.material().source(), &reference);
        assert_eq!(
            sample.material().rheology().temperature().temperature_k(),
            temperature(1.4)
        );
        near(flow.pressure_pa()[cell], pressure(1.4), 4e-5);
    }
    assert_eq!(previous, saved);
}

#[test]
fn g4_adaptive_work_source_boundary_and_cancellation_refusals_are_atomic() {
    let gate = CancelGate::new_clock_free();
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let reference = transport(2.0, 0.4);
    let owner = MaterialFluxSystem::new(
        &system,
        &reference,
        MaterialFluxNumerics {
            linear_tolerance: 1e-12,
            ..MaterialFluxNumerics::default()
        },
    )
    .unwrap();
    let previous = moments(&mesh, &|_| [1.0, 0.0]);
    let saved = previous.clone();
    let force = |time: f64, _: [f64; 2]| [6.0 * time * time, 0.0];
    let traction = |_: f64, _: [f64; 2], _: [f64; 2]| [0.0; 2];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|edge| mesh.is_boundary(*edge))
        .map(|edge| MaterialNavierStokesTimeBoundary::Open {
            edge,
            traction: &traction,
            backflow: OpenBoundaryPolicy::NaturalTraction,
        })
        .collect();
    let mut limited = controls(1e-12);
    limited.maximum_attempts = 1;
    let work = owner.advance_navier_stokes_mixed_rheology_sdirk2(
        0.0,
        0.4,
        &previous,
        &force,
        &boundaries,
        &limited,
        |_, _, _| {
            Ok(MaterialFluxViscosityTangent {
                material: reference,
                derivative_pa_s2: 0.0,
            })
        },
        &gate,
    );
    assert!(matches!(
        work,
        Err(MaterialSdirk2AdaptiveError::WorkLimit {
            attempts: 1,
            last_numerical_failure: None
        })
    ));
    let mut config = controls(1e-6);
    config.initial_step_s = 0.01;
    config.maximum_step_s = 0.01;
    let refusal = MaterialError::State {
        what: "late original temperature support".into(),
    };
    let mut refused_queries = 0;
    let source = owner.advance_navier_stokes_mixed_rheology_sdirk2(
        0.0,
        0.05,
        &previous,
        &force,
        &boundaries,
        &config,
        |time, _, _| {
            if time > 0.025 {
                refused_queries += 1;
                return Err(refusal.clone());
            }
            Ok(MaterialFluxViscosityTangent {
                material: reference,
                derivative_pa_s2: 0.0,
            })
        },
        &gate,
    );
    assert_eq!(refused_queries, 1);
    assert!(
        matches!(source, Err(MaterialSdirk2AdaptiveError::Stage(MaterialSdirk2StageFailure {
        trial: MaterialSdirk2Trial::Coarse, time_s,
        error: FluxSdirk2Error { stage: Some(FluxSdirk2Stage::Endpoint), error: MaterialFluxError::Material(error) },
    })) if error == refusal && (time_s - 0.03).abs() < 1e-15)
    );
    let bad = |time: f64, _: [f64; 2], _: [f64; 2]| {
        if time > 0.025 {
            [f64::NAN, 0.0]
        } else {
            [0.0; 2]
        }
    };
    let bad_boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|edge| mesh.is_boundary(*edge))
        .map(|edge| MaterialNavierStokesTimeBoundary::Open {
            edge,
            traction: &bad,
            backflow: OpenBoundaryPolicy::NaturalTraction,
        })
        .collect();
    let boundary = owner.advance_navier_stokes_mixed_rheology_sdirk2(
        0.0,
        0.05,
        &previous,
        &force,
        &bad_boundaries,
        &config,
        |_, _, _| {
            Ok(MaterialFluxViscosityTangent {
                material: reference,
                derivative_pa_s2: 0.0,
            })
        },
        &gate,
    );
    assert!(matches!(
        boundary,
        Err(MaterialSdirk2AdaptiveError::Stage(
            MaterialSdirk2StageFailure {
                error: FluxSdirk2Error {
                    error: MaterialFluxError::Newton(FluxNewtonError::Boundary(_)),
                    ..
                },
                ..
            }
        ))
    ));
    let retry = owner
        .advance_navier_stokes_mixed_rheology_sdirk2(
            0.0,
            0.05,
            &previous,
            &force,
            &boundaries,
            &config,
            |_, _, _| {
                Ok(MaterialFluxViscosityTangent {
                    material: reference,
                    derivative_pa_s2: 0.0,
                })
            },
            &gate,
        )
        .unwrap();
    let (error, _) = system.velocity_error(
        &retry.second_half.endpoint.flow.step.solution.discrete().x,
        &|_| [1.0 + 0.05_f64.powi(3), 0.0],
    );
    assert!(error < 1e-7);
    let unrepresentable_time = owner.advance_navier_stokes_mixed_rheology_sdirk2(
        1e16,
        1e16 + 2.0,
        &previous,
        &|_, _| panic!("collapsed stage queried force"),
        &boundaries,
        &config,
        |_, _, _| -> Result<MaterialFluxViscosityTangent<FluidTransportState>, MaterialError> {
            panic!("collapsed stage queried source")
        },
        &gate,
    );
    assert!(matches!(
        unrepresentable_time,
        Err(MaterialSdirk2AdaptiveError::StepUnderflow { .. })
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let result = owner.advance_navier_stokes_mixed_rheology_sdirk2(
        0.0,
        0.05,
        &previous,
        &|_, _| panic!("cancelled force"),
        &boundaries,
        &config,
        |_, _, _| -> Result<MaterialFluxViscosityTangent<FluidTransportState>, MaterialError> {
            panic!("cancelled material")
        },
        &cancelled,
    );
    assert!(matches!(
        result,
        Err(MaterialSdirk2AdaptiveError::Cancelled)
    ));
    assert_eq!(previous, saved);
}
