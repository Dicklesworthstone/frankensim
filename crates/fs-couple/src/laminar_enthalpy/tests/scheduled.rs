//! Authored mechanical events and coupled endpoint source/time regressions.

use super::super::adaptive::*;
use super::super::boussinesq::BoussinesqConfig;
use super::super::drive::*;
use super::*;
use fs_flux::OpenBoundaryPolicy;
use std::cell::RefCell;

struct Driven {
    query: LaminarDriveQuery,
    boundaries: Vec<LaminarDrivenBoundary>,
    force: [f64; 2],
    gradient: [f64; 2],
    velocity: [f64; 2],
    inlet: f64,
    temperature: f64,
}
impl LaminarMechanicalDrive for Driven {
    fn boundaries(&self) -> &[LaminarDrivenBoundary] {
        &self.boundaries
    }
    fn force_density(&self, _: [f64; 2]) -> [f64; 2] {
        self.force
    }
    fn velocity(&self, _: usize, _: [f64; 2]) -> [f64; 2] {
        self.velocity
    }
    fn traction(&self, _: usize, point: [f64; 2], normal: [f64; 2]) -> [f64; 2] {
        normal.map(|n| -(17.0 + self.gradient[0] * point[0] + self.gradient[1] * point[1]) * n)
    }
}
fn drive(
    mesh: &TriMesh,
    state: &DepositionSimulation<'_>,
    cells: &[CellId],
    query: LaminarDriveQuery,
    prescribed: bool,
) -> Driven {
    Driven {
        query,
        boundaries: (0..mesh.edges.len())
            .filter(|&edge| mesh.is_boundary(edge))
            .map(|edge| {
                if prescribed {
                    LaminarDrivenBoundary::Velocity { edge }
                } else {
                    LaminarDrivenBoundary::Open {
                        edge,
                        backflow: OpenBoundaryPolicy::NaturalTraction,
                    }
                }
            })
            .collect(),
        force: [0.0; 2],
        gradient: [0.0; 2],
        velocity: [0.0; 2],
        inlet: 300.0,
        temperature: cells
            .iter()
            .map(|&id| state.temperature_k(id).unwrap())
            .sum::<f64>()
            / cells.len() as f64,
    }
}
fn accuracy() -> LaminarAdaptiveControls {
    LaminarAdaptiveControls {
        step: controls(),
        initial_step_s: 0.037,
        minimum_step_s: 1e-5,
        maximum_step_s: 0.037,
        maximum_native_steps: 300,
        maximum_accepted_steps: 100,
        maximum_endpoint_evaluations: 30_000,
        energy_absolute_tolerance_j: 0.05,
        temperature_absolute_tolerance_k: 0.05,
        velocity_absolute_tolerance_m_s: 1e-4,
        transfer_absolute_tolerance_j: 1e-5,
        relative_tolerance: 0.0,
        window_tolerance_j: 1e-6,
    }
}

#[test]
fn opposing_time_dependent_force_and_pressure_work_refine_even_at_constant_velocity_and_temperature()
 {
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
        let velocity = |_| [0.3, 0.0];
        let previous = moments(&mesh, &velocity);
        let target = 0.111;
        let solved = owner
            .advance_scheduled(
                cx,
                &old,
                target,
                LaminarEnthalpyHistory {
                    velocity_m_s: &previous,
                    time_s: 0.0,
                },
                &[],
                LaminarAdaptiveControls {
                    initial_step_s: target,
                    maximum_step_s: target,
                    ..accuracy()
                },
                |_, query, state| {
                    let mut source = drive(&mesh, state, &cells, query, false);
                    source.force = [1.0 + 2.0 * query.time_s, 0.0];
                    source.gradient = source.force;
                    Ok::<_, ()>(source)
                },
                |_, t, _, _| Ok(source(t, 0.4)),
                |_, _, _, _| Ok::<_, ()>(exchange(&mesh, law.h(300.0))),
            )
            .unwrap();
        let report = solved.report();
        assert!(report.work.rejected_comparisons > 0);
        assert!(report.estimates.iter().all(|e| e.energy_error_ratio < 0.01
            && e.temperature_error_ratio < 0.01
            && e.velocity_error_ratio < 0.01));
        assert!(
            report
                .estimates
                .iter()
                .any(|e| e.transfer_error_ratio > 0.2)
        );
        let expected: f64 = report
            .steps
            .iter()
            .map(|s| DEPTH * 0.3 * (1.0 + 2.0 * s.end_s) * (s.end_s - s.start_s))
            .sum();
        let exact = DEPTH * 0.3 * (target + target * target);
        let mechanical = report.mechanical.unwrap();
        assert!(mechanical.body_work_j > exact);
        near(mechanical.body_work_j, expected, 2e-8);
        near(mechanical.imposed_traction_work_j, -expected, 2e-8);
        near(
            mechanical.body_work_j + mechanical.imposed_traction_work_j,
            0.0,
            2e-8,
        );
        near(mechanical.outward_kinetic_energy_j, 0.0, 2e-9);
        near(report.stored_energy_change_j, 0.0, 2e-7);
        near(report.thermal_balance_residual_j, 0.0, 1e-6);
        assert_eq!(report.events_reached, 0);
        let (error, divergence) =
            system.velocity_error(&solved.flow().step.solution.discrete().x, &velocity);
        assert!(error < 2e-8 && divergence < 2e-8);
    });
}

#[test]
fn decimal_velocity_events_reverse_original_streams_keep_left_limits_and_resume_the_same_physics() {
    with_cx(|cx, _| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 1e5 };
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
        let events = [0.037, 0.074, 0.111];
        let speeds = [0.3, -0.2, 0.1];
        let inlets = [320.0, 300.0, 330.0];
        // This fixture isolates exact stops and native stream recurrence;
        // the preceding test independently demands work-driven refinement.
        let event_accuracy = LaminarAdaptiveControls {
            transfer_absolute_tolerance_j: 0.05,
            ..accuracy()
        };
        let previous = moments(&mesh, &|_| [speeds[0], 0.0]);
        let history = LaminarEnthalpyHistory {
            velocity_m_s: &previous,
            time_s: 0.0,
        };
        let queries = RefCell::new(Vec::new());
        let resolve = |_: &Cx<'_>, query: LaminarDriveQuery, state: &DepositionSimulation<'_>| {
            queries.borrow_mut().push(query);
            assert!(query.segment < 3);
            let lower = if query.segment == 0 {
                0.0
            } else {
                events[query.segment - 1]
            };
            assert!(query.time_s > lower && query.time_s <= events[query.segment]);
            let mut source = drive(&mesh, state, &cells, query, true);
            source.velocity = [speeds[query.segment], 0.0];
            source.inlet = inlets[query.segment];
            Ok::<_, ()>(source)
        };
        let thermal_exchange =
            |_: &Cx<'_>,
             _: &DepositionSimulation<'_>,
             _: &MaterialNavierStokesRheologyStep<
                LaminarThermalMaterial<ThermalFluidTransportState>,
            >,
             source: &Driven| Ok::<_, ()>(exchange(&mesh, law.h(source.inlet)));
        let solved = owner
            .advance_scheduled(
                cx,
                &old,
                events[2],
                history,
                &events,
                event_accuracy,
                &resolve,
                |_, t, _, _| Ok(source(t, 0.4)),
                &thermal_exchange,
            )
            .unwrap();
        let report = solved.report();
        assert_eq!(report.events_reached, 3);
        assert_eq!(report.steps.len(), 6);
        for (segment, &time_s) in events.iter().enumerate() {
            assert!(
                queries
                    .borrow()
                    .contains(&LaminarDriveQuery { time_s, segment })
            );
        }
        let mut temperatures = [298.0, 308.0];
        let mut heat = 0.0;
        let mut wall_work = 0.0;
        let mut previous_speed = speeds[0];
        for step in &report.steps {
            let source = &step.report.exchange.drive;
            assert_eq!(source.query, step.report.exchange.query);
            assert_eq!(source.query.time_s, step.end_s);
            let speed = speeds[source.query.segment];
            assert_eq!(source.velocity[0], speed);
            wall_work += RHO * DEPTH * (speed - previous_speed) * speed;
            previous_speed = speed;
            let first = mesh
                .edges
                .iter()
                .find(|edge| edge.tris.1 == usize::MAX && speed * edge.normal[0] < -0.01)
                .unwrap()
                .tris
                .0;
            let second = 1 - first;
            let dt = step.end_s - step.start_s;
            let a = 2.0 * speed.abs() * dt;
            temperatures[first] = (temperatures[first] + a * source.inlet) / (1.0 + a);
            temperatures[second] = (temperatures[second] + a * temperatures[first]) / (1.0 + a);
            heat += RHO * DEPTH * speed.abs() * CP * dt * (source.inlet - temperatures[second]);
            near(
                source.temperature,
                0.5 * (temperatures[0] + temperatures[1]),
                2e-7,
            );
        }
        for cell in 0..2 {
            near(
                solved.thermal().temperature_k(cells[cell]).unwrap(),
                temperatures[cell],
                2e-7,
            );
        }
        near(report.external_stream_enthalpy_j, heat, 3e-7);
        near(report.stored_energy_change_j, heat, 3e-7);
        near(
            report.mechanical.unwrap().velocity_boundary_work_j.unwrap(),
            wall_work,
            3e-8,
        );
        near(wall_work, 0.078, 1e-12);
        let first = owner
            .advance_scheduled(
                cx,
                &old,
                events[0],
                history,
                &events,
                event_accuracy,
                &resolve,
                |_, t, _, _| Ok(source(t, 0.4)),
                &thermal_exchange,
            )
            .unwrap();
        let resumed = owner
            .advance_scheduled(
                cx,
                first.thermal(),
                events[2],
                LaminarEnthalpyHistory {
                    velocity_m_s: first.flow().step.solution.velocity_dofs(),
                    time_s: events[0],
                },
                &events,
                event_accuracy,
                &resolve,
                |_, t, _, _| Ok(source(t, 0.4)),
                &thermal_exchange,
            )
            .unwrap();
        assert_eq!(resumed.report().events_reached, 2);
        assert_eq!(resumed.report().steps[0].report.exchange.query.segment, 1);
        for &cell in &cells {
            near(
                resumed.thermal().temperature_k(cell).unwrap(),
                solved.thermal().temperature_k(cell).unwrap(),
                2e-8,
            );
        }
        near(
            first.report().external_stream_enthalpy_j + resumed.report().external_stream_enthalpy_j,
            heat,
            3e-7,
        );
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(previous, moments(&mesh, &|_| [speeds[0], 0.0]));
    });
}

#[test]
fn scheduled_thermal_actuator_and_original_buoyancy_close_independent_bdf_work_without_heat_duplication()
 {
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
        let events = [0.037, 0.074, 0.111];
        let forces = [0.2, -0.4, 0.1];
        let previous = vec![0.0; system.n_u];
        let solved = owner
            .advance_boussinesq_scheduled(
                cx,
                &old,
                events[2],
                LaminarEnthalpyHistory {
                    velocity_m_s: &previous,
                    time_s: 0.0,
                },
                &events,
                BoussinesqConfig {
                    gravity_m_s2: [0.0, -10.0],
                    max_relative_density_departure: 0.05,
                    max_relative_specific_heat_departure: 0.0,
                },
                LaminarAdaptiveControls {
                    velocity_absolute_tolerance_m_s: 1e-5,
                    transfer_absolute_tolerance_j: 1e-6,
                    ..accuracy()
                },
                |_, query, state| {
                    let mut response = drive(&mesh, state, &cells, query, false);
                    response.force[0] =
                        forces[query.segment] * (1.0 + 0.01 * (response.temperature - 300.0));
                    response.inlet = 300.0 + 10.0 * query.time_s;
                    Ok::<_, ()>(response)
                },
                |_, t, _, _| {
                    ThermalFluidTransportState::try_new(
                        FluidTransportState::try_new(
                            RHO * (1.0 - 0.002 * (t - 300.0)),
                            CP,
                            0.5,
                            0.4,
                        )?,
                        t,
                        None,
                    )
                },
                |_, _, _, drive| {
                    let mut response = exchange(&mesh, law.h(drive.inlet));
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
        let (mut velocity, mut body, mut buoyancy, mut defect) = ([0.0_f64; 2], 0.0, 0.0, 0.0);
        for step in &report.steps {
            let query = step.report.exchange.query;
            let dt = step.end_s - step.start_s;
            let force = [
                forces[query.segment] * (1.0 + 0.1 * step.end_s),
                0.4 * step.end_s,
            ];
            let increment = force.map(|f| dt * f / RHO);
            for axis in 0..2 {
                velocity[axis] += increment[axis];
            }
            body += DEPTH * dt * force[0] * velocity[0];
            buoyancy += DEPTH * dt * force[1] * velocity[1];
            defect += 0.5 * RHO * DEPTH * (increment[0].powi(2) + increment[1].powi(2));
            near(
                step.report.exchange.drive.temperature,
                300.0 + 10.0 * step.end_s,
                2e-7,
            );
        }
        assert_eq!(report.events_reached, 3);
        assert!(report.work.rejected_comparisons > 0);
        let (error, divergence) =
            system.velocity_error(&solved.flow().step.solution.discrete().x, &|_| velocity);
        assert!(error < 2e-8 && divergence < 2e-8, "{error}, {divergence}");
        let mechanical = report.mechanical.unwrap();
        near(mechanical.body_work_j, body, 2e-8);
        near(mechanical.imposed_traction_work_j, 0.0, 2e-8);
        near(mechanical.outward_kinetic_energy_j, 0.0, 2e-9);
        near(report.buoyancy_work_j, buoyancy, 2e-8);
        let kinetic = 0.5 * RHO * DEPTH * (velocity[0].powi(2) + velocity[1].powi(2));
        near(body + buoyancy, kinetic + defect, 2e-10);
        assert!(defect > 0.0);
        near(report.stored_energy_change_j, 24.0 * events[2], 2e-7);
        near(report.other_external_heat_j, 24.0 * events[2], 2e-7);
        near(report.bulk_viscous_heat_j, 0.0, 1e-12);
        near(report.thermal_balance_residual_j, 0.0, 1e-6);
        assert_eq!(old.time_s(), 0.0);
        assert!(previous.iter().all(|&v| v == 0.0));
    });
}

#[test]
fn prescribed_acceleration_work_refines_when_endpoint_velocity_and_temperature_agree() {
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
        let (initial, acceleration, target) = (0.1, 0.5, 0.2);
        let previous = moments(&mesh, &|_| [initial, 0.0]);
        let solved = owner
            .advance_scheduled(
                cx,
                &old,
                target,
                LaminarEnthalpyHistory {
                    velocity_m_s: &previous,
                    time_s: 0.0,
                },
                &[],
                LaminarAdaptiveControls {
                    initial_step_s: target,
                    maximum_step_s: target,
                    transfer_absolute_tolerance_j: 1e-5,
                    ..accuracy()
                },
                |_, query, state| {
                    let mut value = drive(&mesh, state, &cells, query, true);
                    value.velocity = [initial + acceleration * query.time_s, 0.0];
                    Ok::<_, ()>(value)
                },
                |_, t, _, _| Ok(source(t, 0.4)),
                |_, _, _, _| Ok::<_, ()>(exchange(&mesh, law.h(300.0))),
            )
            .unwrap();
        let report = solved.report();
        assert!(report.work.rejected_comparisons > 0);
        // The original imposed velocity gives identical coarse/fine endpoint
        // states. Only the new prescribed-work channel requests refinement.
        assert!(report.estimates.iter().all(|e| e.energy_error_ratio < 0.01
            && e.temperature_error_ratio < 0.01
            && e.velocity_error_ratio < 0.01));
        assert!(
            report
                .estimates
                .iter()
                .any(|e| e.transfer_error_ratio > 0.2)
        );
        let (mut work, mut damping) = (0.0, 0.0);
        for step in &report.steps {
            let speed = initial + acceleration * step.end_s;
            let increment = acceleration * (step.end_s - step.start_s);
            let expected = RHO * DEPTH * increment * speed;
            work += expected;
            damping += 0.5 * RHO * DEPTH * increment * increment;
            near(
                step.report
                    .mechanical
                    .unwrap()
                    .velocity_boundary_work_j
                    .unwrap(),
                expected,
                2e-8,
            );
        }
        let final_speed = initial + acceleration * target;
        let kinetic_gain = 0.5 * RHO * DEPTH * (final_speed * final_speed - initial * initial);
        near(work, kinetic_gain + damping, 1e-12);
        assert!(damping > 0.0);
        let mechanical = report.mechanical.unwrap();
        near(mechanical.velocity_boundary_work_j.unwrap(), work, 3e-8);
        near(mechanical.body_work_j, 0.0, 1e-12);
        near(mechanical.imposed_traction_work_j, 0.0, 1e-12);
        near(mechanical.outward_kinetic_energy_j, 0.0, 1e-10);
        near(mechanical.backflow_resistance_work_j.unwrap(), 0.0, 1e-12);
        near(report.stored_energy_change_j, 0.0, 1e-7);
        near(report.bulk_viscous_heat_j, 0.0, 1e-12);
        near(report.other_external_heat_j, 0.0, 1e-12);
        near(report.thermal_balance_residual_j, 0.0, 1e-6);
        let (error, divergence) = system
            .velocity_error(&solved.flow().step.solution.discrete().x, &|_| {
                [final_speed, 0.0]
            });
        assert!(error < 2e-8 && divergence < 2e-8);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(previous, moments(&mesh, &|_| [initial, 0.0]));
    });
}

#[test]
fn late_mechanical_event_refusal_cancellation_and_shared_budget_publish_no_prefix() {
    with_cx(|cx, gate| {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let law = Caloric { datum: 0.0 };
        let (old, cells) = inventory(&mesh, &law, &[300.0; 2], false);
        let before = old.cells().clone();
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
        let previous = vec![0.0; system.n_u];
        let history = LaminarEnthalpyHistory {
            velocity_m_s: &previous,
            time_s: 0.0,
        };
        let events = [0.037, 0.074, 0.111];
        let error = owner
            .advance_scheduled(
                cx,
                &old,
                events[2],
                history,
                &events,
                LaminarAdaptiveControls {
                    maximum_native_steps: 3,
                    ..accuracy()
                },
                |_, query, state| Ok::<_, ()>(drive(&mesh, state, &cells, query, true)),
                |_, t, _, _| Ok(source(t, 0.4)),
                |_, _, _, _| Ok::<_, ()>(exchange(&mesh, law.h(300.0))),
            )
            .err()
            .unwrap();
        assert!(matches!(error.failure, LaminarAdaptiveFailure::WorkBudget));
        assert_eq!(error.work.native_steps, 3);
        let error = owner
            .advance_scheduled(
                cx,
                &old,
                events[2],
                history,
                &events,
                accuracy(),
                |_, query, state| {
                    if query.segment == 2 {
                        Err("scheduled actuator refused")
                    } else {
                        Ok(drive(&mesh, state, &cells, query, true))
                    }
                },
                |_, t, _, _| Ok(source(t, 0.4)),
                |_, _, _, _| Ok::<_, &str>(exchange(&mesh, law.h(300.0))),
            )
            .err()
            .unwrap();
        assert!(matches!(
            error.failure,
            LaminarAdaptiveFailure::Native {
                error: LaminarEnthalpyError::Drive("scheduled actuator refused"),
                ..
            }
        ));
        assert_eq!(error.work.native_steps, 7);
        let error = owner
            .advance_scheduled(
                cx,
                &old,
                events[2],
                history,
                &events,
                accuracy(),
                |_, query, state| {
                    if query.segment == 2 {
                        gate.request();
                    }
                    Ok::<_, ()>(drive(&mesh, state, &cells, query, true))
                },
                |_, t, _, _| Ok(source(t, 0.4)),
                |_, _, _, _| Ok::<_, ()>(exchange(&mesh, law.h(300.0))),
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
        assert_eq!(error.work.native_steps, 7);
        assert_eq!(old.cells(), &before);
        assert_eq!(old.time_s(), 0.0);
        assert!(previous.iter().all(|&v| v == 0.0));
    });
}
