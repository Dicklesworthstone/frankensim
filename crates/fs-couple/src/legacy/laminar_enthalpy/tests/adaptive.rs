//! Independent finite-volume, acceleration and transfer-quadrature oracles.

use super::super::adaptive::*;
use super::super::boussinesq::BoussinesqConfig;
use super::*;
use fs_flux::OpenBoundaryPolicy;

fn adaptive_controls(tolerance: f64) -> LaminarAdaptiveControls {
    LaminarAdaptiveControls {
        step: controls(),
        initial_step_s: 0.3,
        minimum_step_s: 1e-5,
        maximum_step_s: 0.3,
        maximum_native_steps: 1_000,
        maximum_accepted_steps: 100,
        maximum_endpoint_evaluations: 30_000,
        energy_absolute_tolerance_j: tolerance,
        temperature_absolute_tolerance_k: tolerance,
        velocity_absolute_tolerance_m_s: tolerance,
        transfer_absolute_tolerance_j: tolerance,
        relative_tolerance: 0.0,
        window_tolerance_j: 1e-6,
    }
}

fn input<'a>(
    velocity: &'a [f64],
    boundaries: &'a [MaterialNavierStokesBoundary<'a>],
) -> LaminarEnthalpyInput<'a> {
    LaminarEnthalpyInput {
        history: LaminarEnthalpyHistory {
            velocity_m_s: velocity,
            time_s: 0.0,
        },
        force_density: &zero,
        boundaries,
        slip_heat_routes: &[],
    }
}

#[test]
fn throughflow_refinement_tracks_original_cell_recurrence_and_continuous_limit() {
    with_cx(|cx, _| {
        for datum in [0.0, 1e5] {
            for speed in [-0.7_f64, 0.7] {
                let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
                let system = FluxSystem::new(&mesh);
                let law = Caloric { datum };
                let (old, cells) = inventory(&mesh, &law, &[298.0, 308.0], false);
                let material = reference(0.4);
                let owner = LaminarEnthalpy::new(
                    cx,
                    &system,
                    &material,
                    MaterialFluxNumerics {
                        linear_tolerance: 1e-12,
                        ..MaterialFluxNumerics::default()
                    },
                    &old,
                    map(&cells),
                    DEPTH,
                    2 * mesh.edges.len(),
                )
                .unwrap();
                let velocity = |_| [speed, 0.0];
                let boundaries: Vec<_> = (0..mesh.edges.len())
                    .filter(|&edge| mesh.is_boundary(edge))
                    .map(|edge| MaterialNavierStokesBoundary::Velocity {
                        edge,
                        value: &velocity,
                    })
                    .collect();
                let previous = moments(&mesh, &velocity);
                let original = old.cells().clone();
                let first = mesh
                    .edges
                    .iter()
                    .find(|edge| edge.tris.1 == usize::MAX && speed * edge.normal[0] < -0.1)
                    .unwrap()
                    .tris
                    .0;
                let second = 1 - first;
                let initial = [298.0, 308.0];
                let target = 0.3;
                let lambda = 2.0 * speed.abs();
                let decay = (-lambda * target).exp();
                let exact_first = 330.0 + (initial[first] - 330.0) * decay;
                let exact_second = 330.0
                    + (initial[second] - 330.0 + lambda * target * (initial[first] - 330.0))
                        * decay;
                let mut errors = Vec::new();
                let mut trajectories = Vec::new();
                for tolerance in [0.05, 0.003] {
                    let solved = owner
                        .advance_adaptive(
                            cx,
                            &old,
                            target,
                            input(&previous, &boundaries),
                            adaptive_controls(tolerance),
                            |_, t, _| Ok(source(t, 0.4)),
                            |_, _, _| Ok::<_, ()>(exchange(&mesh, law.h(330.0))),
                        )
                        .unwrap();
                    let report = solved.report();
                    assert!(report.work.rejected_comparisons > 0);
                    assert!((1e-5..=0.3).contains(&report.suggested_step_s));
                    assert_eq!(report.steps.len(), 2 * report.estimates.len());
                    assert_eq!(
                        report.work.native_steps,
                        3 * (report.estimates.len() + report.work.rejected_comparisons),
                    );
                    assert!(
                        report.work.endpoint_evaluations
                            > report
                                .steps
                                .iter()
                                .map(|s| s.report.endpoint_evaluations)
                                .sum::<usize>()
                    );
                    let mut predicted = initial;
                    let mut heat = 0.0;
                    let mut time = 0.0;
                    for step in &report.steps {
                        assert_eq!(step.start_s, time);
                        let dt = step.end_s - step.start_s;
                        let a = lambda * dt;
                        predicted[first] = (predicted[first] + a * 330.0) / (1.0 + a);
                        predicted[second] = (predicted[second] + a * predicted[first]) / (1.0 + a);
                        heat += RHO * DEPTH * speed.abs() * CP * dt * (330.0 - predicted[second]);
                        time = step.end_s;
                    }
                    assert_eq!(time, target);
                    assert_eq!(solved.thermal().time_s(), target);
                    for cell in 0..2 {
                        near(
                            solved.thermal().temperature_k(cells[cell]).unwrap(),
                            predicted[cell],
                            2e-7,
                        );
                        near(
                            solved.flow().samples()[cell]
                                .material()
                                .source()
                                .temperature_k(),
                            predicted[cell],
                            2e-7,
                        );
                    }
                    near(report.external_stream_enthalpy_j, heat, 3e-7);
                    near(report.stored_energy_change_j, heat, 3e-7);
                    near(report.thermal_balance_residual_j, 0.0, 1e-6);
                    near(report.bulk_viscous_heat_j, 0.0, 1e-12);
                    assert_eq!(report.buoyancy_work_j, 0.0);
                    assert!(
                        report
                            .estimates
                            .iter()
                            .all(|e| e.maximum_error_ratio() <= 1.0)
                    );
                    let error = (predicted[first] - exact_first)
                        .abs()
                        .max((predicted[second] - exact_second).abs());
                    errors.push(error);
                    trajectories.push(report.steps.len());
                }
                assert!(errors[1] < 0.6 * errors[0], "{errors:?}");
                assert!(trajectories[1] > trajectories[0]);
                assert_eq!(old.cells(), &original);
                assert_eq!(old.time_s(), 0.0);
                assert_eq!(previous, moments(&mesh, &velocity));
            }
        }
    });
}

#[test]
fn buoyancy_refines_actual_momentum_when_native_uniform_heating_is_already_exact() {
    with_cx(|cx, _| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 1e5 };
        let (old, cells) = inventory(&mesh, &law, &[300.0; 2], false);
        let material = reference(2.0);
        let owner = LaminarEnthalpy::new(
            cx,
            &system,
            &material,
            MaterialFluxNumerics {
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
            &old,
            map(&cells),
            DEPTH,
            2 * mesh.edges.len(),
        )
        .unwrap();
        let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
        // Literal natural traction retains the interior incoming momentum,
        // so uniform acceleration is an exact spatial native solution.
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&edge| mesh.is_boundary(edge))
            .map(|edge| MaterialNavierStokesBoundary::Open {
                edge,
                traction: &traction,
                backflow: OpenBoundaryPolicy::NaturalTraction,
            })
            .collect();
        let previous = vec![0.0; system.n_u];
        let mut accuracy = adaptive_controls(1e-5);
        accuracy.initial_step_s = 0.5;
        accuracy.maximum_step_s = 0.5;
        accuracy.velocity_absolute_tolerance_m_s = 1e-4;
        accuracy.transfer_absolute_tolerance_j = 1.0;
        let config = BoussinesqConfig {
            gravity_m_s2: [0.0, -10.0],
            max_relative_density_departure: 0.05,
            max_relative_specific_heat_departure: 0.0,
        };
        let solved = owner
            .advance_boussinesq_adaptive(
                cx,
                &old,
                1.0,
                input(&previous, &boundaries),
                config,
                accuracy,
                |_, t, _| {
                    ThermalFluidTransportState::try_new(
                        FluidTransportState::try_new(
                            RHO * (1.0 - 0.002 * (t - 300.0)),
                            CP,
                            0.5,
                            2.0,
                        )?,
                        t,
                        None,
                    )
                },
                |_, state, _| {
                    let mut response = exchange(&mesh, law.h(300.0 + 10.0 * state.time_s()));
                    for &cell in &cells {
                        response
                            .heat
                            .bulk
                            .sources_w
                            .insert(cell, old.cells()[&cell].mass_kg * CP * 10.0);
                    }
                    Ok::<_, ()>(response)
                },
            )
            .unwrap();
        let report = solved.report();
        assert!(report.work.rejected_comparisons > 0);
        assert!(
            report
                .estimates
                .iter()
                .any(|e| e.velocity_error_ratio > 0.2)
        );
        assert!(
            report
                .estimates
                .iter()
                .all(|e| e.temperature_error_ratio < 0.01)
        );
        let mut velocity_y = 0.0;
        let mut work_j = 0.0;
        let mut square_spans = 0.0;
        for step in &report.steps {
            let dt = step.end_s - step.start_s;
            // Original rho(T)-rho0 gives a_y(t)=0.2*t. This recurrence
            // independently verifies the accepted velocity history chaining.
            velocity_y += 0.2 * step.end_s * dt;
            square_spans += dt * dt;
            work_j += DEPTH * 0.4 * step.end_s * velocity_y * dt;
            let buoyancy = step.report.buoyancy.as_ref().unwrap();
            near(
                buoyancy.cell_force_density_n_m3[0][1],
                0.4 * step.end_s,
                2e-8,
            );
        }
        let expected = |_| [0.0, velocity_y];
        let (error, divergence) =
            system.velocity_error(&solved.flow().step.solution.discrete().x, &expected);
        assert!(error < 2e-8 && divergence < 2e-8, "{error}, {divergence}");
        near(velocity_y, 0.1 + 0.1 * square_spans, 1e-13);
        near(report.buoyancy_work_j, work_j, 2e-8);
        assert!(work_j > 0.0);
        near(report.stored_energy_change_j, 24.0, 2e-7);
        near(report.other_external_heat_j, 24.0, 2e-7);
        near(report.external_stream_enthalpy_j, 0.0, 2e-7);
        near(report.bulk_viscous_heat_j, 0.0, 1e-12);
        near(report.thermal_balance_residual_j, 0.0, 1e-6);
        for (index, sample) in solved.flow().samples().iter().enumerate() {
            near(
                solved.thermal().temperature_k(cells[index]).unwrap(),
                310.0,
                2e-7,
            );
            near(
                sample
                    .material()
                    .source()
                    .resolved_transport()
                    .density_kg_m3(),
                0.98 * RHO,
                2e-9,
            );
            assert_eq!(sample.transport().density_kg_m3(), RHO);
            assert_eq!(sample.transport().specific_heat_j_kg_k(), CP);
        }
        assert_eq!(old.time_s(), 0.0);
        assert!(previous.iter().all(|&value| value == 0.0));
    });
}

#[test]
fn cancelling_stream_and_cooling_heat_still_require_temporal_transfer_resolution() {
    with_cx(|cx, _| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 1e5 };
        let (old, cells) = inventory(&mesh, &law, &[300.0; 2], false);
        let material = reference(0.4);
        let owner = LaminarEnthalpy::new(
            cx,
            &system,
            &material,
            MaterialFluxNumerics {
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
            &old,
            map(&cells),
            DEPTH,
            2 * mesh.edges.len(),
        )
        .unwrap();
        let speed = 0.4;
        let velocity = |_| [speed, 0.0];
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&edge| mesh.is_boundary(edge))
            .map(|edge| MaterialNavierStokesBoundary::Velocity {
                edge,
                value: &velocity,
            })
            .collect();
        let previous = moments(&mesh, &velocity);
        let first = mesh
            .edges
            .iter()
            .find(|edge| edge.tris.1 == usize::MAX && speed * edge.normal[0] < -0.1)
            .unwrap()
            .tris
            .0;
        let mut accuracy = adaptive_controls(1e-6);
        accuracy.transfer_absolute_tolerance_j = 1e-3;
        let heat_slope = RHO * DEPTH * speed * CP * 10.0;
        let solved = owner
            .advance_adaptive(
                cx,
                &old,
                0.2,
                input(&previous, &boundaries),
                accuracy,
                |_, t, _| Ok(source(t, 0.4)),
                |_, state, _| {
                    let mut response = exchange(&mesh, law.h(300.0 + 10.0 * state.time_s()));
                    response
                        .heat
                        .bulk
                        .sources_w
                        .insert(cells[first], -heat_slope * state.time_s());
                    Ok::<_, ()>(response)
                },
            )
            .unwrap();
        let report = solved.report();
        assert!(report.work.rejected_comparisons > 0);
        assert!(
            report
                .estimates
                .iter()
                .any(|e| e.transfer_error_ratio > 0.2)
        );
        let expected_heat: f64 = report
            .steps
            .iter()
            .map(|s| heat_slope * s.end_s * (s.end_s - s.start_s))
            .sum();
        assert!(expected_heat > 0.5 * heat_slope * 0.2 * 0.2);
        near(report.external_stream_enthalpy_j, expected_heat, 2e-7);
        near(report.other_external_heat_j, -expected_heat, 2e-7);
        near(report.stored_energy_change_j, 0.0, 2e-7);
        near(report.thermal_balance_residual_j, 0.0, 1e-6);
        for &cell in &cells {
            near(solved.thermal().temperature_k(cell).unwrap(), 300.0, 2e-7);
        }
    });
}

#[test]
fn native_work_source_and_late_cancellation_refusals_publish_no_window_prefix() {
    with_cx(|cx, gate| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 0.0 };
        let (old, cells) = inventory(&mesh, &law, &[300.0; 2], false);
        let original = old.cells().clone();
        let material = reference(0.4);
        let owner = LaminarEnthalpy::new(
            cx,
            &system,
            &material,
            MaterialFluxNumerics::default(),
            &old,
            map(&cells),
            DEPTH,
            2 * mesh.edges.len(),
        )
        .unwrap();
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&edge| mesh.is_boundary(edge))
            .map(|edge| MaterialNavierStokesBoundary::Velocity { edge, value: &zero })
            .collect();
        let previous = vec![0.0; system.n_u];
        let mut accuracy = adaptive_controls(1e-4);
        accuracy.initial_step_s = 0.1;
        accuracy.maximum_step_s = 0.1;
        let no_heat = |_: &Cx<'_>,
                       _: &DepositionSimulation<'_>,
                       _: &MaterialNavierStokesRheologyStep<
            LaminarThermalMaterial<ThermalFluidTransportState>,
        >| Ok::<_, ()>(exchange(&mesh, law.h(300.0)));
        let mut budget = accuracy;
        budget.maximum_native_steps = 3;
        let error = owner
            .advance_adaptive(
                cx,
                &old,
                0.3,
                input(&previous, &boundaries),
                budget,
                |_, t, _| Ok(source(t, 0.4)),
                &no_heat,
            )
            .err()
            .unwrap();
        assert!(matches!(error.failure, LaminarAdaptiveFailure::WorkBudget));
        assert_eq!(error.work.native_steps, 3);
        assert!(error.work.endpoint_evaluations > 0);
        budget = accuracy;
        budget.maximum_endpoint_evaluations = 1;
        let error = owner
            .advance_adaptive(
                cx,
                &old,
                0.3,
                input(&previous, &boundaries),
                budget,
                |_, t, _| Ok(source(t, 0.4)),
                &no_heat,
            )
            .err()
            .unwrap();
        assert!(matches!(
            error.failure,
            LaminarAdaptiveFailure::Native {
                error: LaminarEnthalpyError::EndpointBudget,
                endpoint_allowance: 1,
            }
        ));
        assert_eq!(error.work.native_steps, 1);
        assert_eq!(error.work.endpoint_evaluations, 0);
        let calls = Cell::new(0usize);
        let error = owner
            .advance_adaptive(
                cx,
                &old,
                0.3,
                input(&previous, &boundaries),
                accuracy,
                |_, t, _| {
                    calls.set(calls.get() + 1);
                    Ok(source(t, 0.4))
                },
                |_, state, _| {
                    if state.time_s() > 0.15 {
                        Err("authored source refusal")
                    } else {
                        Ok(exchange(&mesh, law.h(300.0)))
                    }
                },
            )
            .err()
            .unwrap();
        assert!(matches!(
            error.failure,
            LaminarAdaptiveFailure::Native {
                error: LaminarEnthalpyError::Exchange("authored source refusal"),
                ..
            }
        ));
        assert_eq!(error.work.native_steps, 4);
        assert!(calls.get() > 0);
        let error = owner
            .advance_adaptive(
                cx,
                &old,
                0.3,
                input(&previous, &boundaries),
                accuracy,
                |_, t, _| Ok(source(t, 0.4)),
                |_, state, _| {
                    if state.time_s() > 0.15 {
                        gate.request();
                    }
                    Ok::<_, ()>(exchange(&mesh, law.h(300.0)))
                },
            )
            .err()
            .unwrap();
        assert!(matches!(
            error.failure,
            LaminarAdaptiveFailure::Native {
                error: LaminarEnthalpyError::Cancelled,
                ..
            }
        ));
        assert_eq!(error.work.native_steps, 4);
        assert_eq!(old.cells(), &original);
        assert_eq!(old.time_s(), 0.0);
        assert!(previous.iter().all(|&value| value == 0.0));
    });
}

#[test]
fn original_inventory_events_and_unresolved_time_floor_refuse_without_accepting_a_prefix() {
    with_cx(|cx, _| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 0.0 };
        let (plain, cells) = inventory(&mesh, &law, &[300.0; 2], false);
        let mut events: Vec<_> = plain
            .cells()
            .values()
            .cloned()
            .map(|cell| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            })
            .collect();
        let mut boundary = plain.cells()[&cells[0]].boundary;
        boundary.source_w = 1.0;
        events.push(TimedEvent {
            time_s: 0.15,
            event: Event::Boundary(cells[0], boundary),
        });
        let old = DepositionSimulation::new(
            events,
            &[(MATERIAL, &law as &dyn ThermalLaw)],
            plain.budget(),
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 1e-10,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
        let material = reference(0.4);
        let owner = LaminarEnthalpy::new(
            cx,
            &system,
            &material,
            MaterialFluxNumerics::default(),
            &old,
            map(&cells),
            DEPTH,
            2 * mesh.edges.len(),
        )
        .unwrap();
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&edge| mesh.is_boundary(edge))
            .map(|edge| MaterialNavierStokesBoundary::Velocity { edge, value: &zero })
            .collect();
        let previous = vec![0.0; system.n_u];
        let error = owner
            .advance_adaptive(
                cx,
                &old,
                0.3,
                input(&previous, &boundaries),
                adaptive_controls(0.01),
                |_, t, _| Ok(source(t, 0.4)),
                |_, _, _| Ok::<_, ()>(exchange(&mesh, law.h(300.0))),
            )
            .err()
            .unwrap();
        assert_eq!(error.work.native_steps, 0);
        assert!(matches!(
            error.failure,
            LaminarAdaptiveFailure::Native {
                error: LaminarEnthalpyError::Thermal(DepositionError::InvalidInput(_)),
                ..
            }
        ));
        // Three decimal caps must hit the exact target without a one-ULP
        // tail demanding an impossible extra coarse/fine comparison.
        let mut decimal = adaptive_controls(0.01);
        decimal.initial_step_s = 0.037;
        decimal.maximum_step_s = 0.037;
        let rounded = owner
            .advance_adaptive(
                cx,
                &old,
                0.111,
                input(&previous, &boundaries),
                decimal,
                |_, t, _| Ok(source(t, 0.4)),
                |_, _, _| Ok::<_, ()>(exchange(&mesh, law.h(300.0))),
            )
            .unwrap();
        assert_eq!(rounded.thermal().time_s(), 0.111);
        assert_eq!(rounded.report().steps.len(), 6);
        assert_eq!(rounded.report().work.native_steps, 9);
        // Exactly targeting the event returns its original LEFT limit.
        let accepted = owner
            .advance_adaptive(
                cx,
                &old,
                0.15,
                input(&previous, &boundaries),
                adaptive_controls(0.01),
                |_, t, _| Ok(source(t, 0.4)),
                |_, _, _| Ok::<_, ()>(exchange(&mesh, law.h(300.0))),
            )
            .unwrap();
        assert_eq!(accepted.thermal().cells()[&cells[0]].boundary.source_w, 0.0);
        let error = owner
            .advance_adaptive(
                cx,
                accepted.thermal(),
                0.2,
                LaminarEnthalpyInput {
                    history: LaminarEnthalpyHistory {
                        velocity_m_s: accepted.flow().step.solution.velocity_dofs(),
                        time_s: 0.15,
                    },
                    ..input(&previous, &boundaries)
                },
                adaptive_controls(0.01),
                |_, t, _| Ok(source(t, 0.4)),
                |_, _, _| Ok::<_, ()>(exchange(&mesh, law.h(300.0))),
            )
            .err()
            .unwrap();
        assert_eq!(error.work.native_steps, 0);
        assert!(matches!(
            error.failure,
            LaminarAdaptiveFailure::Native {
                error: LaminarEnthalpyError::Thermal(DepositionError::InvalidInput(_)),
                ..
            }
        ));
        let mut accuracy = adaptive_controls(1e-5);
        accuracy.initial_step_s = 0.1;
        accuracy.minimum_step_s = 0.1;
        accuracy.maximum_step_s = 0.1;
        let error = owner
            .advance_adaptive(
                cx,
                &old,
                0.1,
                input(&previous, &boundaries),
                accuracy,
                |_, t, _| Ok(source(t, 0.4)),
                |_, state, _| {
                    let mut response = exchange(&mesh, law.h(300.0));
                    for &cell in &cells {
                        response
                            .heat
                            .bulk
                            .sources_w
                            .insert(cell, 10.0 * state.time_s());
                    }
                    Ok::<_, ()>(response)
                },
            )
            .err()
            .unwrap();
        assert!(matches!(
            error.failure,
            LaminarAdaptiveFailure::AccuracyFloor(_)
        ));
        assert_eq!(error.work.native_steps, 3);
        assert_eq!(error.work.rejected_comparisons, 1);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.cells(), plain.cells());
    });
}
