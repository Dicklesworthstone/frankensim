//! Native pipe reconstruction through prescribed and coupled flow owners.

use super::*;
use crate::graph::GraphSolveConfig;
use crate::graph::inertial::InertialGraphControls;
use crate::graph::thermal::transport::transient::coupled::{
    QuasiSteadyGraphControls, QuasiSteadyGraphProblem, QuasiSteadyGraphResponse,
    advance_quasi_steady_graph_distributed,
};
use crate::graph::thermal::transport::transient::inertial::{
    InertialGraphHeat, InertialThermalControls, InertialThermalProblem, ResolvedInertialGraph,
    advance_inertial_graph_distributed,
};
use crate::graph::thermal::transport::transient::nonlinear::{
    ImplicitTransportControls, ImplicitTransportError,
};
use fs_qty::VolumetricFlowRate;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::{BoundaryStepCoefficients, StepCoefficients};
use std::cell::Cell;

// This file is itself a `#[path]` module of `distributed/tests.rs`, so a bare `mod explicit;` looks beside it in `tests/`, not in `tests/reconstruction/`.
#[path = "reconstruction/explicit.rs"]
mod explicit;

fn controls() -> ImplicitTransportControls {
    ImplicitTransportControls {
        solver: NewtonKrylovConfig {
            absolute_tolerance: 2e-11,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        maximum_cells: 64,
        maximum_newton_iterations: 48,
        maximum_endpoint_evaluations: 10_000,
        energy_scale_j: 100.0,
        finite_difference_step: 1e-5,
        cell_tolerance_j: 1e-6,
        total_tolerance_j: 1e-6,
    }
}

fn no_heat() -> BoundaryStepCoefficients<()> {
    StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    }
    .into()
}

fn wave_average(start: f64, end: f64, displacement: f64) -> f64 {
    let k = 2.0 * core::f64::consts::PI;
    let half = 0.5 * k * (end - start);
    330.0 + 20.0 * (k * (0.5 * (start + end) - displacement)).cos() * half.sin() / half
}

fn wave<'a>(law: &'a Caloric, n: usize) -> (DepositionSimulation<'a>, Vec<CellId>) {
    let (mut data, path) = cells(law, n, 1.0, 1.0 / n as f64);
    for cell in &mut data {
        let t = wave_average(cell.geometry.start_m[0], cell.geometry.end_m[0], 0.0);
        cell.specific_h_j_kg = law.datum + law.cp * (t - 300.0);
    }
    let events = data
        .into_iter()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: n + 2,
            max_links: n + 4,
            max_events: n + 2,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-6,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    (thermal, path)
}

fn inlet_rise(n: usize, q: f64, time: f64) -> f64 {
    30.0 + 20.0 * (2.0 * core::f64::consts::PI * (-1.0 / n as f64 - q * time)).cos()
}

#[test]
fn g1_g3_native_limited_pipe_preserves_heat_and_reduces_smooth_transport_error() {
    with_cx(|cx| {
        let law = law(1e6);
        let n = 16;
        let (initial, path) = wave(&law, n);
        let graph = graph(cx, 1.0);
        let base = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(gas().density),
            &NODES,
            &[path.clone()],
            &initial,
            32,
        )
        .unwrap();
        let limited = base
            .clone()
            .with_limited_reconstruction(cx, &graph, &initial, n - 2)
            .unwrap();
        assert_eq!(limited.reconstruction().unwrap().len(), n - 2);
        assert!(base.reconstruction().is_none());
        assert_eq!(limited.branch_volumes_m3(), base.branch_volumes_m3());
        let actual_inertia = limited.inertial_branches()[0].unwrap();
        let original_inertia = base.inertial_branches()[0].unwrap();
        assert_eq!(actual_inertia.length_m, original_inertia.length_m);
        assert_eq!(actual_inertia.flow_area_m2, original_inertia.flow_area_m2);
        let q = 0.2;
        let dt = 0.01;
        let mut errors = Vec::new();
        for inventory in [&base, &limited] {
            let mut thermal = initial.clone();
            for step in 1..=25 {
                let target = step as f64 * dt;
                let binding = inventory
                    .bind(
                        cx,
                        &graph,
                        &[q],
                        &thermal,
                        streams(&law, q, NODES[0], NODES[1], inlet_rise(n, q, target)),
                        32,
                    )
                    .unwrap();
                let (next, report) = binding
                    .advance_implicit_with_heat(
                        cx,
                        &thermal,
                        target,
                        controls(),
                        &vec![0.0; n + 2],
                        |_, _| Ok::<_, &'static str>(no_heat()),
                    )
                    .unwrap();
                assert!(report.maximum_cell_residual_j <= controls().cell_tolerance_j);
                assert!(report.thermal.residual_j.abs() <= controls().total_tolerance_j);
                for &id in next.cells().keys() {
                    let t = next.temperature_k(id).unwrap();
                    assert!(
                        (310.0 - 1e-7..=350.0 + 1e-7).contains(&t),
                        "bounded transported heat: {t}"
                    );
                }
                thermal = next;
            }
            // Exact cell averages of the continuous plug-flow equation. Keep
            // inlet/outlet and the deliberately first-order end cells outside
            // this smooth-interior spatial-accuracy comparison.
            let error: f64 = path[3..n - 3]
                .iter()
                .map(|id| {
                    let geometry = thermal.cells()[id].geometry;
                    (thermal.temperature_k(*id).unwrap()
                        - wave_average(
                            geometry.start_m[0],
                            geometry.end_m[0],
                            q * thermal.time_s(),
                        ))
                    .abs()
                })
                .sum();
            errors.push(error);
        }
        assert!(
            errors[1] < 0.8 * errors[0],
            "limited/upwind smooth-interior errors: {errors:?}"
        );
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[test]
fn g1_native_quasi_steady_and_retained_momentum_share_the_limited_endpoint() {
    with_cx(|cx| {
        let law = law(0.0);
        let n = 5;
        let (initial, path) = wave(&law, n);
        let graph = graph(cx, 1.0);
        let inventory = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(gas().density),
            &NODES,
            &[path],
            &initial,
            16,
        )
        .unwrap()
        .with_limited_reconstruction(cx, &graph, &initial, n - 2)
        .unwrap();
        let q = 0.2;
        let target = 0.05;
        let boundaries = || streams(&law, q, NODES[0], NODES[1], inlet_rise(n, q, target));
        let binding = inventory
            .bind(cx, &graph, &[q], &initial, boundaries(), 16)
            .unwrap();
        let (prescribed, _) = binding
            .advance_implicit_with_heat(cx, &initial, target, controls(), &[0.0; 7], |_, _| {
                Ok::<_, &'static str>(no_heat())
            })
            .unwrap();
        let pressure = [boundary(0, q * q), boundary(1, 0.0)];
        // `GraphSolveConfig` has no `Default` (its budgets are explicit); the sweep and bisection budgets are the shared graph fixture's.
        let hydraulic = GraphSolveConfig {
            absolute_flow_tolerance: VolumetricFlowRate::new(1e-14),
            relative_flow_tolerance: 1e-12,
            ..crate::graph::tests::config()
        };
        let c = controls();
        let qs_controls = QuasiSteadyGraphControls {
            solver: c.solver,
            maximum_cells: c.maximum_cells,
            maximum_newton_iterations: c.maximum_newton_iterations,
            maximum_endpoint_evaluations: c.maximum_endpoint_evaluations,
            energy_scale_j: c.energy_scale_j,
            finite_difference_step: c.finite_difference_step,
            cell_tolerance_j: c.cell_tolerance_j,
            total_tolerance_j: c.total_tolerance_j,
            maximum_neglected_power_w: 1.0,
            maximum_neglected_energy_j: 1.0,
        };
        let (quasi, quasi_report) = advance_quasi_steady_graph_distributed(
            cx,
            &QuasiSteadyGraphProblem {
                topology: &graph,
                node_cells: &NODES,
                reference_density: Density::new(gas().density),
                maximum_connections: 16,
            },
            &inventory,
            &initial,
            target,
            qs_controls,
            &[0.0; 7],
            |cx, _| {
                Ok::<_, &'static str>(QuasiSteadyGraphResponse {
                    hydraulics: graph
                        .solve(&pressure, hydraulic, cx)
                        .map_err(|_| "native hydraulic solve")?,
                    mass_boundaries: boundaries(),
                    heat: no_heat(),
                })
            },
        )
        .unwrap();
        near(quasi_report.hydraulics.branches[0].flow.value(), q, 1e-13);
        let momentum = InertialGraph::try_new(
            cx,
            &graph,
            Density::new(gas().density),
            &[inventory.inertial_branches()[0].unwrap()],
            &pressure,
            &[q],
            0.0,
            8,
        )
        .unwrap();
        let inertial_controls = InertialThermalControls {
            solver: c.solver,
            maximum_cells: c.maximum_cells,
            maximum_newton_iterations: c.maximum_newton_iterations,
            maximum_endpoint_evaluations: c.maximum_endpoint_evaluations,
            energy_scale_j: c.energy_scale_j,
            flow_scale_m3_s: q,
            finite_difference_step: c.finite_difference_step,
            cell_tolerance_j: c.cell_tolerance_j,
            total_thermal_tolerance_j: c.total_tolerance_j,
            impulse_tolerance_pa_s: 1e-8,
            mechanical_tolerance_j: 1e-8,
            average_flow_tolerance_m3_s: 1e-9,
            momentum: InertialGraphControls {
                maximum_drag_step: 2.0,
                maximum_neglected_power_w: 1.0,
                maximum_neglected_energy_j: 1.0,
            },
        };
        let (next_momentum, inertial, report) = advance_inertial_graph_distributed(
            cx,
            &InertialThermalProblem {
                momentum: &momentum,
                node_cells: &NODES,
                maximum_connections: 16,
            },
            &inventory,
            &initial,
            target,
            inertial_controls,
            &[0.0; 7],
            &[q],
            |_, _, _| {
                Ok::<_, &'static str>(ResolvedInertialGraph {
                    graph: graph.clone(),
                    pressure_boundaries: pressure.to_vec(),
                    driving_pressure_pa: vec![0.0],
                    payload: (),
                })
            },
            |_, _, _, _| {
                Ok::<_, &'static str>(InertialGraphHeat {
                    mass_boundaries: boundaries(),
                    heat: no_heat(),
                })
            },
        )
        .unwrap();
        near(report.hydraulics.average_volume_flows_m3_s[0], q, 1e-9);
        for &id in initial.cells().keys() {
            near(
                quasi.temperature_k(id).unwrap(),
                prescribed.temperature_k(id).unwrap(),
                1e-7,
            );
            near(
                inertial.temperature_k(id).unwrap(),
                prescribed.temperature_k(id).unwrap(),
                1e-7,
            );
        }
        assert_eq!(next_momentum.accepted_time_s(), inertial.time_s());
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(momentum.accepted_time_s(), 0.0);
    });
}

#[test]
fn g0_g4_limited_explicit_mode_retains_linear_implicit_refusal_and_failed_source_rollback() {
    with_cx(|cx| {
        let law = law(0.0);
        let n = 5;
        let (initial, path) = wave(&law, n);
        let graph = graph(cx, 1.0);
        let base = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(gas().density),
            &NODES,
            &[path],
            &initial,
            16,
        )
        .unwrap();
        assert!(
            base.clone()
                .with_limited_reconstruction(cx, &graph, &initial, n - 3)
                .is_err()
        );
        let limited = base
            .with_limited_reconstruction(cx, &graph, &initial, n - 2)
            .unwrap();
        for q in [0.2, -0.2, 0.0] {
            let external = if q == 0.0 {
                vec![]
            } else if q > 0.0 {
                streams(&law, q, NODES[0], NODES[1], 40.0)
            } else {
                streams(&law, q, NODES[1], NODES[0], 40.0)
            };
            let binding = limited
                .bind(cx, &graph, &[q], &initial, external, 16)
                .unwrap();
            assert!(binding.endpoint_transport(cx, &initial).is_ok());
            assert!(binding.wrap_observer(()).is_ok());
            let (explicit, energy, transport) = binding
                .advance_to(
                    cx,
                    &initial,
                    0.1,
                    StepConfig {
                        max_dt_s: 0.01,
                        energy_tolerance_j: 1e-6,
                    },
                )
                .unwrap();
            assert_eq!(explicit.time_s(), 0.1);
            near(energy.residual_j, 0.0, 1e-6);
            near(energy.external_heat_j, transport.external_enthalpy_j, 1e-6);
            for (&id, old) in initial.cells() {
                let temperature = explicit.temperature_k(id).unwrap();
                assert!((310.0 - 1e-8..=350.0 + 1e-8).contains(&temperature));
                if q == 0.0 {
                    assert_eq!(explicit.cells()[&id].specific_h_j_kg, old.specific_h_j_kg);
                }
            }
            assert!(
                binding
                    .advance_advection_implicit(cx, &initial, 0.1, 1e-6, 1e-6)
                    .is_err()
            );
            let calls = Cell::new(0);
            let failure = binding.advance_implicit_with_heat(
                cx,
                &initial,
                0.1,
                controls(),
                &[0.0; 7],
                |_, _| {
                    calls.set(calls.get() + 1);
                    Err::<BoundaryStepCoefficients<()>, _>("source refused")
                },
            );
            assert!(matches!(
                failure,
                Err(ImplicitTransportError::Resolve("source refused"))
            ));
            assert_eq!(calls.get(), 1);
        }
        let binding = limited
            .bind(
                cx,
                &graph,
                &[0.2],
                &initial,
                streams(&law, 0.2, NODES[0], NODES[1], 40.0),
                16,
            )
            .unwrap();
        let calls = Cell::new(0);
        let failure = binding.advance_implicit_with_heat(
            cx,
            &initial,
            0.1,
            ImplicitTransportControls {
                maximum_endpoint_evaluations: 1,
                ..controls()
            },
            &[0.0; 7],
            |_, _| {
                calls.set(calls.get() + 1);
                Ok::<_, &'static str>(no_heat())
            },
        );
        assert!(matches!(
            failure,
            Err(ImplicitTransportError::WorkBudget {
                performed: 1,
                maximum: 1
            })
        ));
        assert_eq!(calls.get(), 1);
        assert_eq!(initial.time_s(), 0.0);
    });
}
