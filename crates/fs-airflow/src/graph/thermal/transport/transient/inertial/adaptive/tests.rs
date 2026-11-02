//! Original coupled owners, with an independent continuous-time finite-volume
//! oracle. The oracle tests temporal convergence, not continuum mesh accuracy.
use super::super::advance_inertial_graph_distributed;
use super::*;
use crate::graph::LossGraph;
use crate::graph::inertial::InertialGraphControls;
use crate::graph::porous::{LinearLossElement, LinearLossResistance, PassiveBranch};
use crate::graph::tests::{boundary, with_cx};
use crate::{SourceProvenance, ToleranceBasis};
use fs_qty::Density;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, Event, MaterialId, StepConfig,
    ThermalLaw, TimedEvent,
};
use std::collections::BTreeMap;

const NODES: [CellId; 2] = [CellId(0), CellId(1)];
const ROUTE: [CellId; 6] = [
    CellId(0),
    CellId(3),
    CellId(2),
    CellId(1),
    CellId(5),
    CellId(4),
];

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h / 1000.0;
        if t.is_finite() && (200.0..=500.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.001
    }
}

fn fixture<'m>(
    cx: &Cx<'_>,
    law: &'m Caloric,
) -> (
    LossGraph,
    DepositionSimulation<'m>,
    DistributedGraphInventory,
    InertialGraph,
) {
    let graph = LossGraph::new_mixed(
        2,
        [(0, 1), (1, 0)]
            .into_iter()
            .enumerate()
            .map(|(i, (from, to))| PassiveBranch::Linear {
                from,
                to,
                loss: LinearLossElement::new(
                    format!("reversing-ring-{i}"),
                    LinearLossResistance::new(1.0),
                    0.0,
                    SourceProvenance::new("analytic six-cell ring", "resolution-test-v1"),
                    ToleranceBasis::EngineeringAllowance,
                )
                .unwrap(),
            })
            .collect(),
    )
    .unwrap();
    let geometry = [
        ([0.0, 0.0, 0.0], [0.0, 0.1, 0.0]),
        ([2.0, 0.0, 0.0], [2.0, 0.1, 0.0]),
        ([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]),
        ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
        ([2.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
    ];
    let schedule = geometry
        .into_iter()
        .enumerate()
        .map(|(i, (start_m, end_m))| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: 0,
                material_id: MaterialId(1),
                mass_kg: 1.0,
                specific_h_j_kg: if i == 0 { 60_000.0 } else { 0.0 },
                geometry: CellGeometry {
                    start_m,
                    end_m,
                    volume_m3: 1.0,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    let thermal = DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 6,
            max_links: 8,
            max_events: 6,
            max_substeps: 256,
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
    let layout = DistributedGraphInventory::new(
        cx,
        &graph,
        Density::new(1.0),
        &NODES,
        &[vec![CellId(3), CellId(2)], vec![CellId(5), CellId(4)]],
        &thermal,
        8,
    )
    .unwrap();
    let sections: Vec<_> = layout
        .inertial_branches()
        .iter()
        .map(|v| v.unwrap())
        .collect();
    let momentum = InertialGraph::try_new(
        cx,
        &graph,
        layout.reference_density(),
        &sections,
        &[boundary(0, 0.0)],
        &[1.0; 2],
        0.0,
        8,
    )
    .unwrap();
    (graph, thermal, layout, momentum)
}

fn step_controls() -> InertialThermalControls {
    InertialThermalControls {
        solver: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        maximum_cells: 8,
        maximum_newton_iterations: 40,
        maximum_endpoint_evaluations: 2048,
        energy_scale_j: 1000.0,
        flow_scale_m3_s: 1.0,
        finite_difference_step: 1e-5,
        cell_tolerance_j: 1e-6,
        total_thermal_tolerance_j: 1e-6,
        impulse_tolerance_pa_s: 1e-8,
        mechanical_tolerance_j: 1e-8,
        average_flow_tolerance_m3_s: 1e-9,
        momentum: InertialGraphControls {
            maximum_drag_step: 2.0,
            maximum_neglected_power_w: 100.0,
            maximum_neglected_energy_j: 100.0,
        },
    }
}

fn limits() -> InertialResolutionControls {
    InertialResolutionControls {
        maximum_interval_s: 1.0,
        minimum_interval_s: 1e-7,
        maximum_native_steps: 2048,
        maximum_accepted_steps: 512,
        maximum_endpoint_evaluations: 200_000,
        absolute_cell_error_j: 100.0,
        relative_cell_error: 0.0,
        absolute_flow_error_m3_s: 1e-3,
        relative_flow_error: 0.0,
        maximum_reversal_step_s: 0.05,
        maximum_cancelled_volume_m3: 1e-4,
    }
}

fn response(graph: &LossGraph, time: f64) -> ResolvedInertialGraph<f64> {
    ResolvedInertialGraph {
        graph: graph.clone(),
        pressure_boundaries: vec![boundary(0, 0.0)],
        driving_pressure_pa: vec![-4.0; 2],
        payload: time,
    }
}
fn heat<T>(payload: T) -> InertialGraphHeat<T> {
    InertialGraphHeat {
        mass_boundaries: vec![],
        heat: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload,
        }
        .into(),
    }
}
fn unchanged(thermal: &DepositionSimulation<'_>, momentum: &InertialGraph) {
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(momentum.accepted_time_s(), 0.0);
    assert_eq!(momentum.volume_flows_m3_s(), vec![1.0; 2]);
    for (&id, cell) in thermal.cells() {
        assert_eq!(
            cell.specific_h_j_kg,
            if id == CellId(0) { 60_000.0 } else { 0.0 }
        );
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}

/// Exact continuous-time solution of the SAME six-tank upwind spatial model:
/// q(t)=5 exp(-t/2)-4. Positive and negative periodic shift operators commute;
/// their exponentials are Poisson kernels. This has no timestep recurrence.
fn continuous_ring_temperature() -> [f64; 6] {
    let zero = 2.0 * 1.25_f64.ln();
    let integral = |t: f64| 10.0 * (1.0 - (-0.5 * t).exp()) - 4.0 * t;
    let positive = integral(zero);
    let negative = positive - integral(1.0);
    let mut result = [300.0; 6];
    let mut forward_weight = 1.0;
    for forward in 0..=20usize {
        let mut backward_weight = 1.0;
        for backward in 0..=20usize {
            let index = (forward + 24 - backward) % 6;
            result[index] += 60.0 * (-positive - negative).exp() * forward_weight * backward_weight;
            backward_weight *= negative / (backward + 1) as f64;
        }
        forward_weight *= positive / (forward + 1) as f64;
    }
    result
}

#[test]
fn g1_reversal_moves_heat_both_ways_when_one_native_average_is_zero_and_refines_in_time() {
    with_cx(|cx| {
        let law = Caloric;
        let (graph, old, layout, momentum) = fixture(cx, &law);
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &NODES,
            maximum_connections: 8,
        };
        let (_, coarse, single) = advance_inertial_graph_distributed(
            cx,
            &problem,
            &layout,
            &old,
            1.0,
            step_controls(),
            &[0.0; 6],
            &[1.0; 2],
            |_, state, _| Ok::<_, &str>(response(&graph, state.time_s())),
            |_, _, _, _| Ok(heat(())),
        )
        .unwrap();
        // I=2, R=1, drive=-4, dt=1: Qnew=-1 and Qaverage=0.
        for q in &single.hydraulics.average_volume_flows_m3_s {
            near(*q, 0.0, 1e-10);
        }
        for (&id, cell) in old.cells() {
            near(
                coarse.cells()[&id].specific_h_j_kg,
                cell.specific_h_j_kg,
                1e-6,
            );
        }

        let calls = Cell::new(0usize);
        let run = |settings| {
            advance_inertial_graph_resolved(
                cx,
                &problem,
                Some(&layout),
                &old,
                1.0,
                step_controls(),
                settings,
                |_, candidate, _| {
                    calls.set(calls.get() + 1);
                    Ok::<_, &str>(response(&graph, candidate.time_s()))
                },
                |_, candidate, &source_time, hydraulic| {
                    assert_eq!(candidate.time_s(), source_time);
                    Ok(heat((
                        candidate.time_s(),
                        hydraulic.average_volume_flows_m3_s[0],
                    )))
                },
            )
        };
        let (_, loose, loose_report) = run(limits()).unwrap();
        calls.set(0);
        let mut tight = limits();
        tight.absolute_cell_error_j = 5.0;
        tight.absolute_flow_error_m3_s = 1e-5;
        tight.maximum_reversal_step_s = 0.005;
        tight.maximum_cancelled_volume_m3 = 1e-6;
        let (fine_momentum, fine, report) = run(tight).unwrap();
        let reference = continuous_ring_temperature();
        near(reference.iter().map(|t| t - 300.0).sum(), 60.0, 1e-11);
        let error = |state: &DepositionSimulation<'_>| {
            ROUTE
                .into_iter()
                .zip(reference)
                .map(|(id, expected)| (state.temperature_k(id).unwrap() - expected).abs())
                .fold(0.0f64, f64::max)
        };
        assert!(error(&fine) < 0.5 * error(&loose));
        assert!(error(&fine) < 0.1);
        assert!(report.steps.len() > loose_report.steps.len());
        // Both sides of the pulse warm, although the coarse signed displacement is zero.
        assert!(fine.temperature_k(ROUTE[1]).unwrap() > 305.0);
        assert!(fine.temperature_k(ROUTE[5]).unwrap() > 305.0);
        for q in fine_momentum.volume_flows_m3_s() {
            near(q, 5.0 * (-0.5_f64).exp() - 4.0, 2e-5);
        }
        assert_eq!(fine.time_s(), 1.0);
        assert_eq!(fine_momentum.accepted_time_s(), 1.0);
        assert_eq!(fine_momentum.branch_count(), 2);
        assert_eq!(fine_momentum.node_count(), 2);
        for (&id, original) in old.cells() {
            let retained = &fine.cells()[&id];
            assert_eq!(retained.mass_kg, original.mass_kg);
            assert_eq!(retained.material_id, original.material_id);
            assert_eq!(retained.geometry, original.geometry);
        }
        assert_eq!(report.work.endpoint_evaluations, calls.get());
        let accepted_calls: usize = report.steps.iter().map(|r| r.endpoint_evaluations).sum();
        assert!(report.work.endpoint_evaluations > accepted_calls);
        assert_eq!(
            report.work.native_steps,
            3 * (report.estimates.len() + report.work.rejected_intervals)
        );
        assert_eq!(report.steps.len(), 2 * report.estimates.len());
        let mut clock = 0.0;
        for step in &report.steps {
            clock += step.hydraulics.duration_s;
            assert_eq!(step.sources, clock);
            assert_eq!(step.heat.0, clock);
            assert_eq!(step.heat.1, step.hydraulics.average_volume_flows_m3_s[0]);
            assert!(step.maximum_cell_residual_j <= step_controls().cell_tolerance_j);
            assert!(
                step.hydraulics.maximum_impulse_residual_pa_s
                    <= step_controls().impulse_tolerance_pa_s
            );
            assert!(
                step.hydraulics.mechanical_balance_residual_j.abs()
                    <= step_controls().mechanical_tolerance_j
            );
            near(
                step.thermal.residual_j,
                0.0,
                step_controls().total_thermal_tolerance_j,
            );
            near(step.thermal.stored_energy_change_j, 0.0, 1e-6);
            near(step.thermal.external_heat_j, 0.0, 1e-12);
        }
        assert_eq!(clock, 1.0);
        assert!(report.estimates.iter().all(|e| e.error_ratio <= 1.0
            && e.maximum_reversal_step_s <= tight.maximum_reversal_step_s
            && e.maximum_cancelled_volume_m3 <= tight.maximum_cancelled_volume_m3));
        unchanged(&old, &momentum);
    });
}

#[test]
fn g0_work_and_resolution_limits_drop_all_local_states_including_accepted_prefixes() {
    with_cx(|cx| {
        let law = Caloric;
        let (graph, old, layout, momentum) = fixture(cx, &law);
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &NODES,
            maximum_connections: 8,
        };
        let calls = Cell::new(0usize);
        let run = |settings| {
            advance_inertial_graph_resolved(
                cx,
                &problem,
                Some(&layout),
                &old,
                1.0,
                step_controls(),
                settings,
                |_, state, _| {
                    calls.set(calls.get() + 1);
                    Ok::<_, &str>(response(&graph, state.time_s()))
                },
                |_, _, _, _| Ok(heat(())),
            )
        };
        let mut native = limits();
        native.maximum_native_steps = 2;
        match run(native) {
            Err(InertialResolutionError::Budget {
                kind: InertialResolutionBudget::NativeSteps,
                performed: 2,
                work,
                ..
            }) => {
                assert_eq!(work.native_steps, 2);
                assert_eq!(work.endpoint_evaluations, calls.get());
            }
            Err(other) => panic!("unexpected native limit refusal: {other:?}"),
            Ok(_) => panic!("native limit accepted"),
        }
        unchanged(&old, &momentum);
        calls.set(0);
        let mut endpoint = limits();
        endpoint.maximum_endpoint_evaluations = 1;
        match run(endpoint) {
            Err(InertialResolutionError::Step {
                source:
                    InertialThermalError::WorkBudget {
                        performed: 1,
                        maximum: 1,
                    },
                work,
                ..
            }) => {
                assert_eq!(work.native_steps, 1);
                assert_eq!(work.endpoint_evaluations, 1);
                assert_eq!(calls.get(), 1);
            }
            Err(other) => panic!("unexpected endpoint limit refusal: {other:?}"),
            Ok(_) => panic!("endpoint limit accepted"),
        }
        unchanged(&old, &momentum);
        calls.set(0);
        let mut accepted = limits();
        accepted.maximum_interval_s = 0.01;
        accepted.maximum_accepted_steps = 2;
        match run(accepted) {
            Err(InertialResolutionError::Budget {
                kind: InertialResolutionBudget::AcceptedSteps,
                performed: 2,
                work,
                ..
            }) => {
                // One pair was locally accepted, but neither original owner advances.
                assert_eq!(work.native_steps, 3);
                assert_eq!(work.rejected_intervals, 0);
                assert_eq!(work.endpoint_evaluations, calls.get());
            }
            Err(other) => panic!("unexpected accepted limit refusal: {other:?}"),
            Ok(_) => panic!("accepted-step limit ignored"),
        }
        unchanged(&old, &momentum);
        let mut floor = limits();
        floor.minimum_interval_s = 1.0;
        match run(floor) {
            Err(InertialResolutionError::Underresolved { estimate, work }) => {
                assert!(estimate.error_ratio > 1.0);
                assert_eq!(work.native_steps, 3);
                assert_eq!(work.rejected_intervals, 1);
            }
            Err(other) => panic!("unexpected resolution-floor refusal: {other:?}"),
            Ok(_) => panic!("underresolved interval accepted"),
        }
        unchanged(&old, &momentum);
    });
}

#[test]
fn g0_resolve_and_exchange_refusals_are_terminal_even_in_discarded_comparison_steps() {
    with_cx(|cx| {
        let law = Caloric;
        let (graph, old, layout, momentum) = fixture(cx, &law);
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &NODES,
            maximum_connections: 8,
        };
        let calls = Cell::new(0usize);
        let result = advance_inertial_graph_resolved(
            cx,
            &problem,
            Some(&layout),
            &old,
            1.0,
            step_controls(),
            limits(),
            |_, _, _| {
                calls.set(calls.get() + 1);
                Err::<ResolvedInertialGraph<()>, _>("unsupported material law")
            },
            |_, _, _, _| Ok(heat(())),
        );
        match result {
            Err(InertialResolutionError::Step {
                source: InertialThermalError::Resolve("unsupported material law"),
                work,
                ..
            }) => {
                assert_eq!(work.native_steps, 1);
                assert_eq!(work.endpoint_evaluations, 1);
                assert_eq!(work.rejected_intervals, 0);
            }
            Err(other) => panic!("unexpected source refusal: {other:?}"),
            Ok(_) => panic!("unsupported source accepted"),
        }
        assert_eq!(calls.get(), 1);
        unchanged(&old, &momentum);
        calls.set(0);
        let result = advance_inertial_graph_resolved(
            cx,
            &problem,
            Some(&layout),
            &old,
            1.0,
            step_controls(),
            limits(),
            |_, state, _| {
                calls.set(calls.get() + 1);
                Ok::<_, &str>(response(&graph, state.time_s()))
            },
            |_, _, _, _| Err::<InertialGraphHeat<()>, _>("unsupported flow regime"),
        );
        assert!(matches!(
            result,
            Err(InertialResolutionError::Step {
                source: InertialThermalError::Exchange("unsupported flow regime"),
                work: InertialResolutionWork {
                    native_steps: 1,
                    endpoint_evaluations: 1,
                    rejected_intervals: 0
                },
                ..
            })
        ));
        assert_eq!(calls.get(), 1);
        unchanged(&old, &momentum);
    });
}

#[test]
fn g0_cancellation_during_original_source_refresh_is_transactional() {
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
        let law = Caloric;
        let (graph, old, layout, momentum) = fixture(&cx, &law);
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &NODES,
            maximum_connections: 8,
        };
        let calls = Cell::new(0usize);
        let exchanges = Cell::new(0usize);
        let result = advance_inertial_graph_resolved(
            &cx,
            &problem,
            Some(&layout),
            &old,
            1.0,
            step_controls(),
            limits(),
            |_, state, _| {
                calls.set(calls.get() + 1);
                gate.request();
                Ok::<_, &str>(response(&graph, state.time_s()))
            },
            |_, _, _, _| {
                exchanges.set(exchanges.get() + 1);
                Ok(heat(()))
            },
        );
        assert!(matches!(
            result,
            Err(InertialResolutionError::Step {
                source: InertialThermalError::Thermal(DepositionError::Cancelled),
                work: InertialResolutionWork {
                    native_steps: 1,
                    endpoint_evaluations: 1,
                    rejected_intervals: 0
                },
                ..
            })
        ));
        assert_eq!(calls.get(), 1);
        assert_eq!(exchanges.get(), 0);
        unchanged(&old, &momentum);
    });
}
