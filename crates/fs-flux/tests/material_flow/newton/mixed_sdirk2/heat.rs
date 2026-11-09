//! G1/G3/G4: original-stage physical heat, refinement and accepted history.

use super::*;
use fs_exec::CancelGate;
use fs_flux::material::adaptive::{
    MaterialNavierStokesTimeBoundary, MaterialSdirk2AdaptiveControls, MaterialSdirk2AdaptiveError,
};

#[test]
fn g1_mixed_sdirk2_bulk_and_slip_heat_have_second_order_time_quadrature() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let density = 3.0;
    let mu = 0.6;
    let beta = 1.2;
    let reference = transport(density, mu);
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
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let dt = 1.0 / steps as f64;
        let mut previous = moments(&mesh, &|p| [p[1], 0.0]);
        let mut total_bulk = 0.0;
        let mut total_slip = 0.0;
        for step in 0..steps {
            let start = step as f64 * dt;
            let times = [
                start + FluxSdirk2Stage::First.time_fraction() * dt,
                start + dt,
            ];
            let reservoir = [0, 1].map(|i| move |p: [f64; 2]| [(1.0 + times[i]) * p[1], 0.0]);
            let traction = [0, 1].map(|i| {
                move |_: [f64; 2], n: [f64; 2]| {
                    let tau = mu * (1.0 + times[i]);
                    [tau * n[1], tau * n[0]]
                }
            });
            let wall = [0, 1].map(|i| {
                move |p: [f64; 2]| {
                    let sign = if p[1] < 0.5 { -1.0 } else { 1.0 };
                    [(1.0 + times[i]) * (p[1] + sign * mu / beta), 0.0]
                }
            });
            let boundaries: [Vec<_>; 2] = std::array::from_fn(|i| {
                mesh.edges
                    .iter()
                    .enumerate()
                    .filter(|(edge, _)| mesh.is_boundary(*edge))
                    .map(|(edge, geometry)| {
                        if geometry.normal[1].abs() > 0.5 {
                            MaterialNavierStokesBoundary::NavierSlip {
                                edge,
                                value: &wall[i],
                                tangential_resistance_pa_s_per_m: beta,
                            }
                        } else if geometry.normal[0] < -0.5 {
                            MaterialNavierStokesBoundary::OpenReservoir {
                                edge,
                                traction: &traction[i],
                                velocity: &reservoir[i],
                            }
                        } else {
                            MaterialNavierStokesBoundary::Open {
                                edge,
                                traction: &traction[i],
                                backflow: OpenBoundaryPolicy::NaturalTraction,
                            }
                        }
                    })
                    .collect()
            });
            let result = owner
                .sdirk2_navier_stokes_mixed_rheology_newton_step(
                    &|_, p| [density * p[1], 0.0],
                    &boundaries[0],
                    &boundaries[1],
                    &previous,
                    dt,
                    FluxNewtonControls {
                        residual_tolerance: 1e-11,
                        ..FluxNewtonControls::default()
                    },
                    |_, _, retained| {
                        assert_eq!(retained, &reference);
                        Ok(MaterialFluxViscosityTangent {
                            material: reference,
                            derivative_pa_s2: 0.0,
                        })
                    },
                )
                .unwrap();
            let heat = &result.dissipated_heat;
            assert_eq!(heat.duration_s(), dt);
            // Shear a=1+t is represented exactly in space and by SDIRK2
            // momentum. Bulk power mu*a² and total slip 2*mu²*a²/beta
            // are both 0.6*(1+t)² W/m, so only heat quadrature has time error.
            for (cell, &energy) in heat.cell_bulk_energy_j_per_m().iter().enumerate() {
                near(energy, mesh.areas[cell] * heat.bulk_energy_j_per_m(), 2e-8);
            }
            for entry in heat.boundary_slip_energy() {
                assert_eq!(entry.component, 0);
                let edge = &mesh.edges[entry.edge];
                let a = mesh.verts[edge.verts.0];
                let b = mesh.verts[edge.verts.1];
                let length = (b[0] - a[0]).hypot(b[1] - a[1]);
                if edge.normal[1].abs() > 0.5 {
                    near(
                        entry.energy_j_per_m,
                        0.5 * length * heat.slip_energy_j_per_m(),
                        2e-8,
                    );
                } else {
                    assert_eq!(entry.energy_j_per_m, 0.0);
                }
            }
            total_bulk += heat.bulk_energy_j_per_m();
            total_slip += heat.slip_energy_j_per_m();
            for stage in [&result.first_stage, &result.endpoint] {
                for sample in stage.flow.samples() {
                    assert_eq!(sample.material(), &reference);
                }
            }
            previous
                .copy_from_slice(&result.endpoint.flow.step.solution.discrete().x[..system.n_u]);
        }
        let exact = 1.4;
        let error = (total_bulk - exact).abs();
        near(total_bulk, total_slip, 2e-7);
        let expected_error = mu * (1.5 * core::f64::consts::FRAC_1_SQRT_2 - 1.0) * dt * dt / 3.0;
        near(error, expected_error, 2e-7);
        errors.push(error);
    }
    for pair in errors.windows(2) {
        assert!((pair[0] / pair[1] - 4.0).abs() < 0.03, "{errors:?}");
    }
}

#[test]
fn g4_adaptive_heat_uses_all_accepted_fine_time_and_discards_rejected_trials() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 2, 2));
    let system = FluxSystem::new(&mesh);
    let density = 3.0;
    let mu = 0.6;
    let beta = 1.2;
    let reference = transport(density, mu);
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
    let start = 2.0;
    let end = 2.4;
    // Cubic translation triggers rejected coarse comparisons while steady
    // unit shear supplies nonzero heat. Exact bulk and total slip powers are
    // both 0.6 W/m; refined temporal velocity error perturbs them only slightly.
    let translation = |t: f64| 1.0 + (t - start).powi(3);
    let velocity = |t, p: [f64; 2]| [translation(t) + p[1], 0.0];
    let wall = |t, p: [f64; 2]| {
        let sign = if p[1] < 0.5 { -1.0 } else { 1.0 };
        [translation(t) + p[1] + sign * mu / beta, 0.0]
    };
    let traction = |_: f64, _: [f64; 2], n: [f64; 2]| [mu * n[1], mu * n[0]];
    let resistance = |_: f64| beta;
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
    let previous = moments(&mesh, &|p| velocity(start, p));
    let saved = previous.clone();
    let controls = MaterialSdirk2AdaptiveControls {
        initial_step_s: 0.4,
        minimum_step_s: 1e-8,
        maximum_step_s: 0.4,
        velocity_absolute_tolerance_m_s: 1e-7,
        velocity_relative_tolerance: 0.0,
        maximum_attempts: 500,
        newton: FluxNewtonControls {
            residual_tolerance: 1e-11,
            ..FluxNewtonControls::default()
        },
        ..MaterialSdirk2AdaptiveControls::default()
    };
    let force = |time: f64, _: [f64; 2]| [3.0 * density * (time - start).powi(2), 0.0];
    let gate = CancelGate::new_clock_free();
    let result = owner
        .advance_navier_stokes_mixed_rheology_sdirk2(
            start,
            end,
            &previous,
            &force,
            &boundaries,
            &controls,
            |_, _, retained| {
                assert_eq!(retained, &reference);
                Ok(MaterialFluxViscosityTangent {
                    material: reference,
                    derivative_pa_s2: 0.0,
                })
            },
            &gate,
        )
        .unwrap();
    assert!(result.report.rejected_comparisons > 0);
    assert!(result.report.accepted_comparisons > 1);
    let heat = &result.dissipated_heat;
    assert_eq!(heat.duration_s(), end - start);
    near(heat.bulk_energy_j_per_m(), mu * (end - start), 2e-4);
    near(heat.slip_energy_j_per_m(), mu * (end - start), 2e-4);
    let last_bulk = result.first_half.dissipated_heat.bulk_energy_j_per_m()
        + result.second_half.dissipated_heat.bulk_energy_j_per_m();
    assert!(heat.bulk_energy_j_per_m() > last_bulk + 0.01);
    let sum_cells: f64 = heat.cell_bulk_energy_j_per_m().iter().sum();
    let sum_edges: f64 = heat
        .boundary_slip_energy()
        .iter()
        .map(|edge| edge.energy_j_per_m)
        .sum();
    near(sum_cells, heat.bulk_energy_j_per_m(), 1e-14);
    near(sum_edges, heat.slip_energy_j_per_m(), 1e-14);
    assert_eq!(previous, saved);

    // A late source refusal must publish neither earlier accepted heat nor a
    // partial flow. Retrying the same immutable history starts a fresh integral.
    let refusal = MaterialError::State {
        what: "late heat source refusal".into(),
    };
    let mut late_controls = controls.clone();
    late_controls.initial_step_s = 0.05;
    late_controls.maximum_step_s = 0.05;
    let failed = owner.advance_navier_stokes_mixed_rheology_sdirk2(
        start,
        end,
        &previous,
        &force,
        &boundaries,
        &late_controls,
        |time, _, _| {
            if time > 2.31 {
                return Err(refusal.clone());
            }
            Ok(MaterialFluxViscosityTangent {
                material: reference,
                derivative_pa_s2: 0.0,
            })
        },
        &gate,
    );
    assert!(
        matches!(failed, Err(MaterialSdirk2AdaptiveError::Stage(failure))
        if matches!(&failure.error.error, MaterialFluxError::Material(error) if error == &refusal))
    );
    let retry = owner
        .advance_navier_stokes_mixed_rheology_sdirk2(
            start,
            end,
            &previous,
            &force,
            &boundaries,
            &controls,
            |_, _, _| {
                Ok(MaterialFluxViscosityTangent {
                    material: reference,
                    derivative_pa_s2: 0.0,
                })
            },
            &gate,
        )
        .unwrap();
    near(
        retry.dissipated_heat.bulk_energy_j_per_m(),
        heat.bulk_energy_j_per_m(),
        1e-13,
    );
    near(
        retry.dissipated_heat.slip_energy_j_per_m(),
        heat.slip_energy_j_per_m(),
        1e-13,
    );
    assert_eq!(previous, saved);
}
