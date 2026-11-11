//! Error-controlled original convection, including latent heat and pump ramps.

use super::*;
use crate::graph::thermal::transport::transient::nonlinear::{
    ImplicitTransportControls, ImplicitTransportError,
    adaptive::{ImplicitTransportWindowControls, ImplicitTransportWindowError},
};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::advection::schedule::AdvectionScheduleKnot;
use fs_thermal::enthalpy::backward_euler::{BoundaryStepCoefficients, StepCoefficients};
use std::cell::Cell;

fn native() -> ImplicitTransportControls {
    ImplicitTransportControls {
        solver: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        maximum_cells: 16,
        maximum_newton_iterations: 32,
        maximum_endpoint_evaluations: 256,
        energy_scale_j: 100.0,
        finite_difference_step: 1e-5,
        cell_tolerance_j: 1e-8,
        total_tolerance_j: 1e-8,
    }
}
fn controls(tolerance: f64) -> ImplicitTransportWindowControls {
    ImplicitTransportWindowControls {
        maximum_interval_s: 0.6,
        minimum_interval_s: 1e-7,
        maximum_native_steps: 20_000,
        maximum_accepted_steps: 4096,
        maximum_source_evaluations: 200_000,
        absolute_specific_energy_error_j_kg: tolerance,
        absolute_temperature_error_k: 1000.0,
        relative_state_error: 0.0,
        energy_tolerance_j: 1e-6,
    }
}
fn heat<T>(payload: T) -> BoundaryStepCoefficients<T> {
    StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload,
    }
    .into()
}

#[test]
fn scheduled_residence_response_refines_and_preserves_actual_native_intervals() {
    with_cx(|cx| {
        let law = Law {
            datum: 1e6,
            plateau: false,
        };
        let old = thermal(&law, [1.0; 2], [0.0; 2], false);
        let schedule = AdvectionSchedule::new(
            [(0.0, 0.0), (0.2, 0.0), (0.4, 1.0), (0.6, 1.0)]
                .into_iter()
                .map(|(time_s, flow_multiplier)| AdvectionScheduleKnot {
                    time_s,
                    flow_multiplier,
                    inlet_enthalpies_j_kg: [(20, law.datum + 100.0)].into(),
                })
                .collect(),
            4,
            1,
            || Ok(()),
        )
        .unwrap();
        let transport = binding(cx, &open_flow(cx, true), &old, streams(law.datum))
            .with_schedule(cx, &old, schedule)
            .unwrap();
        let expected = [
            100.0 * (1.0 - (-0.3_f64).exp()),
            100.0 * (1.0 - 1.3 * (-0.3_f64).exp()),
        ];
        let mut errors = Vec::new();
        for tolerance in [0.02, 0.001] {
            let calls = Cell::new(0);
            let (next, report) = transport
                .advance_implicit_adaptive_with_heat(
                    cx,
                    &old,
                    0.6,
                    native(),
                    controls(tolerance),
                    |_, initial, candidate| {
                        calls.set(calls.get() + 1);
                        assert!(candidate.time_s() > initial.time_s());
                        Ok::<_, &'static str>(heat((initial.time_s(), candidate.time_s())))
                    },
                )
                .unwrap();
            assert_eq!(next.time_s(), 0.6);
            assert_eq!(calls.get(), report.work.source_evaluations);
            assert!(report.work.native_steps > report.steps.len());
            assert!(report.work.rejected_intervals > 0);
            for step in &report.steps {
                assert_eq!(step.report.heat, (step.start_s, step.end_s));
                for knot in [0.2, 0.4] {
                    assert!(!(step.start_s < knot && knot < step.end_s));
                }
            }
            assert!(report.estimates.iter().all(|e| e.error_ratio <= 1.0));
            assert!(report.estimates.iter().any(|e| e.end_s == 0.2));
            assert!(report.estimates.iter().any(|e| e.end_s == 0.4));
            let actual = [
                next.cells()[&CellId(0)].specific_h_j_kg - law.datum,
                next.cells()[&CellId(1)].specific_h_j_kg - law.datum,
            ];
            errors.push((actual[0] - expected[0]).abs() + (actual[1] - expected[1]).abs());
            near(
                report.thermal.stored_energy_change_j,
                actual[0] + actual[1],
                1e-7,
            );
            near(
                report.thermal.external_heat_j,
                report
                    .steps
                    .iter()
                    .map(|s| s.report.thermal.external_heat_j)
                    .sum(),
                1e-7,
            );
            near(report.thermal.residual_j, 0.0, 1e-6);
            assert_eq!(report.thermal.substeps, report.steps.len());
            assert!(report.sum_absolute_native_energy_defect_j <= 1e-6);
        }
        assert!(errors[1] < 0.4 * errors[0], "errors={errors:?}");
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.cells()[&CellId(0)].specific_h_j_kg, law.datum);
    });
}

#[test]
fn latent_plateau_and_opposing_cell_changes_cannot_hide_time_error() {
    with_cx(|cx| {
        for datum in [0.0, 1e6] {
            let law = Law {
                datum,
                plateau: true,
            };
            let old = thermal(&law, [1.0; 2], [30.0, 70.0], false);
            let transport = binding(cx, &loop_flow(cx, false), &old, vec![]);
            let (next, report) = transport
                .advance_implicit_adaptive_with_heat(
                    cx,
                    &old,
                    0.6,
                    native(),
                    controls(0.001),
                    |_, _, _| Ok::<_, &'static str>(heat(())),
                )
                .unwrap();
            assert!(report.work.rejected_intervals > 0);
            assert!(
                report
                    .estimates
                    .iter()
                    .all(|e| e.maximum_temperature_difference_k == 0.0)
            );
            let contrast = 20.0 * (-1.2_f64).exp();
            near(
                next.cells()[&CellId(0)].specific_h_j_kg - datum,
                50.0 - contrast,
                0.08,
            );
            near(
                next.cells()[&CellId(1)].specific_h_j_kg - datum,
                50.0 + contrast,
                0.08,
            );
            assert_eq!(next.temperature_k(CellId(0)).unwrap(), 320.0);
            assert_eq!(next.temperature_k(CellId(1)).unwrap(), 320.0);
            near(report.thermal.stored_energy_change_j, 0.0, 1e-7);
            near(report.thermal.external_heat_j, 0.0, 1e-7);
        }
    });
}

#[test]
fn accepted_heat_is_exactly_the_original_two_interval_source_transactions() {
    with_cx(|cx| {
        let law = Law {
            datum: 0.0,
            plateau: false,
        };
        let old = thermal(&law, [1.0; 2], [0.0; 2], false);
        let transport = binding(cx, &open_flow(cx, false), &old, streams(0.0));
        let resolve =
            |_: &Cx<'_>, initial: &DepositionSimulation<'_>, next: &DepositionSimulation<'_>| {
                // Exact interval mean of authored P(t)=10+2t W. A callback that
                // reused the full window duration or origin would supply wrong heat.
                let power = 10.0 + initial.time_s() + next.time_s();
                let mut coefficients = heat((initial.time_s(), next.time_s(), power));
                coefficients.bulk.sources_w.insert(CellId(1), power);
                Ok::<_, &'static str>(coefficients)
            };
        let (middle, first) = transport
            .advance_implicit_with_heat(cx, &old, 0.3, native(), &[0.0; 2], |cx, next| {
                resolve(cx, &old, next)
            })
            .unwrap();
        let (expected, second) = transport
            .advance_implicit_with_heat(cx, &middle, 0.6, native(), &[0.0; 2], |cx, next| {
                resolve(cx, &middle, next)
            })
            .unwrap();
        let (actual, window) = transport
            .advance_implicit_adaptive_with_heat(cx, &old, 0.6, native(), controls(1e6), resolve)
            .unwrap();
        assert_eq!(actual.cells(), expected.cells());
        assert_eq!(window.steps.len(), 2);
        assert_eq!(window.work.native_steps, 3);
        assert_eq!(window.steps[0].report.heat, first.heat);
        assert_eq!(window.steps[1].report.heat, second.heat);
        near(
            window.thermal.external_heat_j,
            first.thermal.external_heat_j + second.thermal.external_heat_j,
            1e-9,
        );
        near(
            window.thermal.stored_energy_change_j,
            first.thermal.stored_energy_change_j + second.thermal.stored_energy_change_j,
            1e-9,
        );
    });
}

#[test]
fn nonadvecting_wall_storage_controls_the_common_window() {
    with_cx(|cx| {
        let law = Law {
            datum: 0.0,
            plateau: false,
        };
        let base = thermal(&law, [1.0; 2], [0.0; 2], false);
        let mut cells: Vec<_> = base.cells().values().cloned().collect();
        let mut wall = cells[0].clone();
        wall.id = CellId(9);
        wall.segment_id = 9;
        wall.body_id = 9;
        wall.specific_h_j_kg = 100.0;
        wall.boundary.ambient_k = 300.0;
        wall.boundary.htc_w_m2_k = 1.0;
        cells.push(wall);
        let old = DepositionSimulation::new(
            cells
                .into_iter()
                .map(|cell| TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(cell),
                })
                .collect(),
            &[(MaterialId(1), &law as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 3,
                max_links: 0,
                max_events: 3,
                max_substeps: 4096,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-8,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
        let graph = LossGraph::new(2, vec![edge("stopped", 0, 1, 1.0)]).unwrap();
        let transport = GraphEnthalpyTransport::from_volume_flows(
            cx,
            &graph,
            &[0.0],
            Density::new(1.0),
            &[CellId(0), CellId(1)],
            &old,
            vec![],
            8,
        )
        .unwrap();
        let mut c = controls(1e6);
        c.absolute_temperature_error_k = 0.001;
        let (next, report) = transport
            .advance_implicit_adaptive_with_heat(cx, &old, 0.6, native(), c, |_, _, _| {
                Ok::<_, &'static str>(heat(()))
            })
            .unwrap();
        assert!(report.work.rejected_intervals > 0);
        assert_eq!(next.cells()[&CellId(0)], old.cells()[&CellId(0)]);
        assert_eq!(next.cells()[&CellId(1)], old.cells()[&CellId(1)]);
        near(
            next.temperature_k(CellId(9)).unwrap(),
            300.0 + 100.0 * (-0.6_f64).exp(),
            0.15,
        );
        near(
            report.thermal.external_heat_j,
            next.cells()[&CellId(9)].specific_h_j_kg - 100.0,
            1e-6,
        );
    });
}

#[test]
fn source_budget_and_resolution_refusals_leave_inputs_unchanged() {
    with_cx(|cx| {
        let law = Law {
            datum: 0.0,
            plateau: false,
        };
        let old = thermal(&law, [1.0; 2], [0.0; 2], false);
        let transport = binding(cx, &open_flow(cx, false), &old, streams(0.0));
        let refused = transport.advance_implicit_adaptive_with_heat(
            cx,
            &old,
            0.6,
            native(),
            controls(0.01),
            |_, _, _| Err::<BoundaryStepCoefficients<()>, _>("original film source refused"),
        );
        assert!(matches!(
            refused,
            Err(ImplicitTransportWindowError::Step {
                source: ImplicitTransportError::Resolve("original film source refused"),
                ..
            })
        ));
        let mut c = controls(0.01);
        c.maximum_source_evaluations = 1;
        assert!(matches!(
            transport.advance_implicit_adaptive_with_heat(
                cx,
                &old,
                0.6,
                native(),
                c,
                |_, _, _| Ok::<_, &'static str>(heat(())),
            ),
            Err(ImplicitTransportWindowError::Step {
                source: ImplicitTransportError::WorkBudget { .. },
                ..
            })
        ));
        c = controls(1e-12);
        c.minimum_interval_s = 0.6;
        assert!(matches!(
            transport.advance_implicit_adaptive_with_heat(
                cx,
                &old,
                0.6,
                native(),
                c,
                |_, _, _| Ok::<_, &'static str>(heat(())),
            ),
            Err(ImplicitTransportWindowError::Underresolved { .. })
        ));
        c = controls(1e6);
        c.maximum_native_steps = 2;
        assert!(matches!(
            transport.advance_implicit_adaptive_with_heat(
                cx,
                &old,
                0.6,
                native(),
                c,
                |_, _, _| Ok::<_, &'static str>(heat(())),
            ),
            Err(ImplicitTransportWindowError::Budget {
                resource: "native step attempts",
                ..
            })
        ));
        let repeated = || {
            transport
                .advance_implicit_adaptive_with_heat(
                    cx,
                    &old,
                    0.6,
                    native(),
                    controls(0.01),
                    |_, _, _| Ok::<_, &'static str>(heat(())),
                )
                .unwrap()
        };
        let a = repeated();
        let b = repeated();
        assert_eq!(a.0.cells(), b.0.cells());
        assert_eq!(a.1.work, b.1.work);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.cells()[&CellId(0)].specific_h_j_kg, 0.0);
        assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 0.0);
    });
}

#[test]
fn native_events_and_accepted_storage_budget_remain_original_owner_limits() {
    with_cx(|cx| {
        let law = Law {
            datum: 0.0,
            plateau: false,
        };
        let base = thermal(&law, [1.0; 2], [0.0; 2], false);
        let build = |event: bool, limit| {
            let mut events: Vec<_> = base
                .cells()
                .values()
                .cloned()
                .map(|cell| TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(cell),
                })
                .collect();
            if event {
                events.push(TimedEvent {
                    time_s: 0.2,
                    event: Event::Boundary(CellId(0), base.cells()[&CellId(0)].boundary),
                });
            }
            DepositionSimulation::new(
                events,
                &[(MaterialId(1), &law as &dyn ThermalLaw)],
                DepositionBudget {
                    max_cells: 2,
                    // The pipe passes one internal transfer between its two cells, and every
                    // transfer counts against max_links (fs-thermal); no assertion here is about
                    // the link budget.
                    max_links: 1,
                    max_events: 3,
                    max_substeps: limit,
                },
            )
            .unwrap()
            .advance_to(
                0.0,
                StepConfig {
                    max_dt_s: 1.0,
                    energy_tolerance_j: 1e-8,
                },
                || Ok(()),
            )
            .unwrap()
            .0
        };
        let old = build(true, 4096);
        let transport = binding(cx, &open_flow(cx, false), &old, streams(0.0));
        let calls = Cell::new(0);
        let failure = transport.advance_implicit_adaptive_with_heat(
            cx,
            &old,
            0.6,
            native(),
            controls(0.01),
            |_, _, _| {
                calls.set(calls.get() + 1);
                Ok::<_, &'static str>(heat(()))
            },
        );
        assert!(matches!(
            failure,
            Err(ImplicitTransportWindowError::Invalid(_))
        ));
        assert_eq!(calls.get(), 0);
        let (left, _) = transport
            .advance_implicit_adaptive_with_heat(
                cx,
                &old,
                0.2,
                native(),
                controls(1e6),
                |_, _, _| Ok::<_, &'static str>(heat(())),
            )
            .unwrap();
        assert_eq!(left.next_event_time_s(), Some(0.2));
        assert_eq!(old.next_event_time_s(), Some(0.2));
        let bounded = build(false, 1);
        let failure = transport.advance_implicit_adaptive_with_heat(
            cx,
            &bounded,
            0.6,
            native(),
            controls(1e6),
            |_, _, _| {
                calls.set(calls.get() + 1);
                Ok::<_, &'static str>(heat(()))
            },
        );
        assert!(matches!(
            failure,
            Err(ImplicitTransportWindowError::Budget {
                resource: "accepted native thermal steps",
                ..
            })
        ));
        assert_eq!(calls.get(), 0);
    });
}

#[test]
fn cancellation_during_original_source_query_has_no_partial_window() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 711,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let law = Law {
            datum: 0.0,
            plateau: false,
        };
        let old = thermal(&law, [1.0; 2], [0.0; 2], false);
        let transport = binding(&cx, &open_flow(&cx, false), &old, streams(0.0));
        let calls = Cell::new(0);
        let result = transport.advance_implicit_adaptive_with_heat(
            &cx,
            &old,
            0.6,
            native(),
            controls(0.01),
            |_, _, _| {
                calls.set(calls.get() + 1);
                gate.request();
                Ok::<_, &'static str>(heat(()))
            },
        );
        assert!(matches!(
            result,
            Err(ImplicitTransportWindowError::Step {
                source: ImplicitTransportError::Thermal(DepositionError::Cancelled),
                ..
            })
        ));
        assert_eq!(calls.get(), 1);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.cells()[&CellId(0)].specific_h_j_kg, 0.0);
    });
}
