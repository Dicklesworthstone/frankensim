//! Independent finite-network roots exercise the actual momentum and heat owners.
mod affine;
use super::*;
use crate::graph::inertial::InertialBranch;
use crate::graph::porous::{LinearLossElement, LinearLossResistance, PassiveBranch};
use crate::graph::tests::{boundary, with_cx};
use crate::{SourceProvenance, ToleranceBasis};
use fs_qty::Density;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, MaterialId, StepConfig,
    ThermalLaw, TimedEvent,
};
use std::collections::BTreeMap;

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
fn inventory<'a>(
    law: &'a Caloric,
    temperatures: &[f64],
    powers: &[f64],
) -> DepositionSimulation<'a> {
    let schedule = temperatures
        .iter()
        .zip(powers)
        .enumerate()
        .map(|(i, (&t, &power))| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: 0,
                material_id: MaterialId(1),
                mass_kg: 1.0,
                specific_h_j_kg: 1000.0 * (t - 300.0),
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 400.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: power,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: temperatures.len(),
            max_links: 16,
            max_events: 16,
            max_substeps: 100,
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
}
fn network(nodes: usize, ends: &[(usize, usize)], resistances: &[f64]) -> LossGraph {
    let branches = ends
        .iter()
        .zip(resistances)
        .enumerate()
        .map(|(i, (&(from, to), &r))| PassiveBranch::Linear {
            from,
            to,
            loss: LinearLossElement::new(
                format!("branch-{i}"),
                LinearLossResistance::new(r),
                0.0,
                SourceProvenance::new(
                    "authored analytic network fixture",
                    "inertial-thermal-test-v1",
                ),
                ToleranceBasis::EngineeringAllowance,
            )
            .unwrap(),
        })
        .collect();
    LossGraph::new_mixed(nodes, branches).unwrap()
}
fn owner(cx: &Cx<'_>, graph: &LossGraph, inertances: &[f64], initial: &[f64]) -> InertialGraph {
    let sections: Vec<_> = inertances
        .iter()
        .map(|&length| InertialBranch {
            length_m: length,
            flow_area_m2: 1.0,
        })
        .collect();
    InertialGraph::try_new(
        cx,
        graph,
        Density::new(1.0),
        &sections,
        &[boundary(0, 0.0)],
        initial,
        0.0,
        16,
    )
    .unwrap()
}
fn controls() -> InertialThermalControls {
    InertialThermalControls {
        solver: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        maximum_cells: 8,
        maximum_newton_iterations: 40,
        maximum_endpoint_evaluations: 1024,
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
fn response(graph: &LossGraph, drives: &[f64]) -> ResolvedInertialGraph<()> {
    ResolvedInertialGraph {
        graph: graph.clone(),
        pressure_boundaries: vec![boundary(0, 0.0)],
        driving_pressure_pa: drives.to_vec(),
        payload: (),
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
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn g0_final_actual_flow_refresh_rejects_a_sharp_supported_drag_law() {
    with_cx(|cx| {
        let graph = network(2, &[(0, 1), (1, 0)], &[1.0; 2]);
        let momentum = owner(cx, &graph, &[1.0; 2], &[0.0; 2]);
        let law = Caloric;
        let old = inventory(&law, &[300.0; 2], &[0.0; 2]);
        let cells = [CellId(0), CellId(1)];
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &cells,
            maximum_connections: 2,
        };
        let guess = 1.0 - 5e-12;
        let width = 1.0 - guess;
        let calls = Cell::new(0);
        let result = advance_inertial_graph(
            cx,
            &problem,
            &old,
            1.0,
            controls(),
            &[0.0; 2],
            &[guess; 2],
            |_, _, average| {
                calls.set(calls.get() + 1);
                // A continuous, supported source law; no regime/domain refusal
                // exposes the stale coefficient. Both drag states are admissible.
                let fraction = ((average[0] - guess) / width).clamp(0.0, 1.0);
                let r = 1.0 + 0.5 * fraction;
                Ok::<_, String>(ResolvedInertialGraph {
                    graph: network(2, &[(0, 1), (1, 0)], &[r; 2]),
                    pressure_boundaries: vec![boundary(0, 0.0)],
                    driving_pressure_pa: vec![3.0; 2],
                    payload: r,
                })
            },
            |_, _, &r, _| {
                assert!((1.0..=1.5).contains(&r));
                Ok(heat(r))
            },
        );
        // Held R=1 gives Qaverage=1, Qnew=2. The rate mismatch alone is
        // already below the solver and physical tolerances. Refreshed R=1.5
        // gives a 0.5 Pa s impulse defect per edge and a 1 J work defect on
        // that SAME trajectory; stepping again would conceal this diagnosis.
        match result {
            Err(InertialThermalError::Momentum(InertialGraphError::MomentumResidual {
                impulse_pa_s,
                energy_j,
            })) => {
                assert!(impulse_pa_s > 0.49);
                assert!(energy_j > 0.99);
            }
            Err(error) => panic!("unexpected refusal: {error:?}"),
            Ok(_) => panic!("stale candidate-flow law was accepted"),
        }
        assert_eq!(calls.get(), 3);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(momentum.accepted_time_s(), 0.0);
        assert!(old.cells().values().all(|cell| cell.specific_h_j_kg == 0.0));
        assert_eq!(momentum.volume_flows_m3_s(), vec![0.0; 2]);
    });
}

#[test]
fn g0_final_refresh_retains_actual_flow_sources_and_charges_the_budget() {
    with_cx(|cx| {
        let graph = network(2, &[(0, 1), (1, 0)], &[1.0; 2]);
        let momentum = owner(cx, &graph, &[1.0; 2], &[0.0; 2]);
        let law = Caloric;
        let old = inventory(&law, &[300.0; 2], &[0.0; 2]);
        let cells = [CellId(0), CellId(1)];
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &cells,
            maximum_connections: 2,
        };
        let calls = Cell::new(0);
        let run = |limits| {
            advance_inertial_graph(
                cx,
                &problem,
                &old,
                1.0,
                limits,
                &[0.0; 2],
                &[1.0 - 5e-12; 2],
                |_, _, average| {
                    calls.set(calls.get() + 1);
                    Ok::<_, String>(ResolvedInertialGraph {
                        graph: graph.clone(),
                        pressure_boundaries: vec![boundary(0, 0.0)],
                        driving_pressure_pa: vec![3.0; 2],
                        payload: average[0],
                    })
                },
                |_, _, &source_rate, hydraulic| {
                    Ok(heat((source_rate, hydraulic.average_volume_flows_m3_s[0])))
                },
            )
        };
        let mut limited = controls();
        limited.maximum_endpoint_evaluations = 2;
        assert!(matches!(
            run(limited),
            Err(InertialThermalError::WorkBudget {
                performed: 2,
                maximum: 2
            })
        ));
        assert_eq!(calls.get(), 2);
        calls.set(0);
        let (next_momentum, next_thermal, report) = run(controls()).unwrap();
        let actual = report.hydraulics.average_volume_flows_m3_s[0];
        assert_eq!(report.sources, actual);
        assert_eq!(report.heat, (actual, actual));
        assert_eq!(report.endpoint_evaluations, 3);
        assert_eq!(calls.get(), 3);
        assert_eq!(report.nonlinear.iterations, 0);
        assert_eq!(next_momentum.accepted_time_s(), 1.0);
        assert_eq!(next_thermal.time_s(), 1.0);
        assert_eq!(momentum.accepted_time_s(), 0.0);
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn g1_branched_startup_matches_independent_momentum_and_enthalpy_roots() {
    with_cx(|cx| {
        let graph = network(2, &[(0, 1), (0, 1), (1, 0)], &[1.0, 2.0, 1.5]);
        let momentum = owner(cx, &graph, &[2.0, 4.0, 3.0], &[0.0; 3]);
        let law = Caloric;
        let old = inventory(&law, &[330.0, 300.0], &[0.0; 2]);
        let cells = [CellId(0), CellId(1)];
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &cells,
            maximum_connections: 3,
        };
        let (next_momentum, next, report) = advance_inertial_graph(
            cx,
            &problem,
            &old,
            1.0,
            controls(),
            &[0.0; 2],
            &[0.0; 3],
            |_, _, _| Ok::<_, String>(response(&graph, &[1.0, 0.0, 0.0])),
            |_, _, _, hydraulic| Ok(heat(hydraulic.average_volume_flows_m3_s.clone())),
        )
        .unwrap();
        // (2 L/dt + R) Qavg = head + nodal pressure difference at zero initial Q.
        // Eliminate the one free-node pressure using Q0 + Q1 = Q2.
        let p = (1.0 / 5.0) / (1.0 / 5.0 + 1.0 / 10.0 + 1.0 / 7.5);
        let expected = [(1.0 - p) / 5.0, -p / 10.0, p / 7.5];
        for (i, q) in expected.into_iter().enumerate() {
            near(report.hydraulics.average_volume_flows_m3_s[i], q, 1e-10);
            near(report.hydraulics.final_flows_m3_s[i], 2.0 * q, 1e-10);
            near(report.heat[i], q, 1e-10);
        }
        // Both counterflow branches return the same total mass as branch zero.
        let difference = 30.0 / (1.0 + 2.0 * expected[0]);
        near(
            next.temperature_k(CellId(0)).unwrap(),
            315.0 + 0.5 * difference,
            1e-7,
        );
        near(
            next.temperature_k(CellId(1)).unwrap(),
            315.0 - 0.5 * difference,
            1e-7,
        );
        near(report.thermal.stored_energy_change_j, 0.0, 1e-6);
        near(report.thermal.external_heat_j, 0.0, 1e-12);
        assert!(
            report.hydraulics.maximum_impulse_residual_pa_s <= controls().impulse_tolerance_pa_s
        );
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(momentum.accepted_time_s(), 0.0);
        assert_eq!(next.time_s(), 1.0);
        assert_eq!(next_momentum.accepted_time_s(), 1.0);
    });
}

#[test]
fn g1_reversal_transports_heat_with_average_flow_not_endpoint_sign() {
    with_cx(|cx| {
        let graph = network(3, &[(0, 1), (1, 2), (2, 0)], &[1.0; 3]);
        let momentum = owner(cx, &graph, &[1.0; 3], &[1.0; 3]);
        let law = Caloric;
        let old = inventory(&law, &[330.0, 300.0, 300.0], &[0.0; 3]);
        let cells = [CellId(0), CellId(1), CellId(2)];
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &cells,
            maximum_connections: 3,
        };
        let (_, next, report) = advance_inertial_graph(
            cx,
            &problem,
            &old,
            1.0,
            controls(),
            &[0.0; 3],
            &[1.0; 3],
            |_, _, _| Ok::<_, String>(response(&graph, &[-1.5; 3])),
            |_, _, _, _| Ok(heat(())),
        )
        .unwrap();
        // Each branch: Qnew=(0.5*Qold-1.5)/1.5=-2/3; Qavg=1/6 > 0.
        for i in 0..3 {
            near(report.hydraulics.final_flows_m3_s[i], -2.0 / 3.0, 1e-10);
            near(
                report.hydraulics.average_volume_flows_m3_s[i],
                1.0 / 6.0,
                1e-10,
            );
        }
        let ratio: f64 = 1.0 / 7.0;
        let excess = (30.0 / (1.0 + 1.0 / 6.0)) / (1.0 - ratio.powi(3));
        for i in 0..3 {
            near(
                next.temperature_k(CellId(i)).unwrap(),
                300.0 + excess * ratio.powi(i as i32),
                1e-7,
            );
        }
        assert!(next.temperature_k(CellId(1)).unwrap() > next.temperature_k(CellId(2)).unwrap());
        near(report.thermal.stored_energy_change_j, 0.0, 1e-6);
    });
}

#[test]
fn g1_temperature_dependent_drag_closes_on_the_accepted_enthalpy_field() {
    with_cx(|cx| {
        let graph = network(2, &[(0, 1), (1, 0)], &[1.0; 2]);
        let momentum = owner(cx, &graph, &[1.0; 2], &[0.0; 2]);
        let law = Caloric;
        let old = inventory(&law, &[330.0, 300.0], &[100.0, -100.0]);
        let cells = [CellId(0), CellId(1)];
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &cells,
            maximum_connections: 2,
        };
        let (_, next, report) = advance_inertial_graph(
            cx,
            &problem,
            &old,
            1.0,
            controls(),
            &[0.0; 2],
            &[0.0; 2],
            |_, candidate, _| {
                let r = 1.0
                    + 0.01
                        * (candidate
                            .temperature_k(CellId(0))
                            .map_err(|e| e.to_string())?
                            - 300.0);
                let resolved = network(2, &[(0, 1), (1, 0)], &[r; 2]);
                Ok::<_, String>(ResolvedInertialGraph {
                    graph: resolved,
                    pressure_boundaries: vec![boundary(0, 0.0)],
                    driving_pressure_pa: vec![1.0; 2],
                    payload: r,
                })
            },
            |_, _, &r, hydraulic| {
                // Test fixture's explicit constitutive support at actual old/new/average rates.
                if !(1.0..=2.0).contains(&r)
                    || hydraulic.speed_envelope_m3_s.iter().any(|q| *q > 1.0)
                {
                    return Err("fixture source support".to_string());
                }
                Ok(heat(r))
            },
        )
        .unwrap();
        let (mut lo, mut hi) = (0.0, 0.5);
        for _ in 0..80 {
            let q = 0.5 * (lo + hi);
            let t = 315.0 + 15.1 / (1.0 + 2.0 * q);
            if q * (3.0 + 0.01 * (t - 300.0)) > 1.0 {
                hi = q;
            } else {
                lo = q;
            }
        }
        let q = 0.5 * (lo + hi);
        let t0 = 315.0 + 15.1 / (1.0 + 2.0 * q);
        near(report.hydraulics.average_volume_flows_m3_s[0], q, 1e-9);
        near(next.temperature_k(CellId(0)).unwrap(), t0, 1e-7);
        near(report.sources, 1.0 + 0.01 * (t0 - 300.0), 1e-9);
        near(report.heat, report.sources, 1e-12);
        near(report.thermal.residual_j, 0.0, 1e-6);
    });
}

#[test]
fn g0_rest_preserves_heat_storage_and_failures_publish_neither_clock() {
    with_cx(|cx| {
        let graph = network(2, &[(0, 1), (1, 0)], &[1.0; 2]);
        let momentum = owner(cx, &graph, &[1.0; 2], &[0.0; 2]);
        let law = Caloric;
        let old = inventory(&law, &[300.0; 2], &[100.0, 0.0]);
        let cells = [CellId(0), CellId(1)];
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &cells,
            maximum_connections: 2,
        };
        let (_, next, report) = advance_inertial_graph(
            cx,
            &problem,
            &old,
            1.0,
            controls(),
            &[0.0; 2],
            &[0.0; 2],
            |_, _, _| Ok::<_, String>(response(&graph, &[0.0; 2])),
            |_, _, _, _| Ok(heat(())),
        )
        .unwrap();
        near(next.temperature_k(CellId(0)).unwrap(), 300.1, 1e-9);
        assert!(report.hydraulics.final_flows_m3_s.iter().all(|q| *q == 0.0));
        near(report.thermal.external_heat_j, 100.0, 1e-10);
        let mut limited = controls();
        limited.maximum_endpoint_evaluations = 1;
        assert!(matches!(
            advance_inertial_graph(
                cx,
                &problem,
                &old,
                1.0,
                limited,
                &[0.0; 2],
                &[0.0; 2],
                |_, _, _| Ok::<_, String>(response(&graph, &[0.0; 2])),
                |_, _, _, _| Ok(heat(())),
            ),
            Err(InertialThermalError::WorkBudget { .. })
        ));
        assert!(matches!(
            advance_inertial_graph(
                cx,
                &problem,
                &old,
                1.0,
                controls(),
                &[0.0; 2],
                &[0.0; 2],
                |_, _, _| Ok::<_, &str>(response(&graph, &[0.0; 2])),
                |_, _, _, _| Err::<InertialGraphHeat<()>, _>("actual-flow source refusal"),
            ),
            Err(InertialThermalError::Exchange("actual-flow source refusal"))
        ));
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(momentum.accepted_time_s(), 0.0);
        assert!(old.cells().values().all(|c| c.specific_h_j_kg == 0.0));
    });
}

#[test]
fn g0_direct_flow_binding_requires_real_mass_reservoirs_at_pressure_nodes() {
    with_cx(|cx| {
        let graph = network(2, &[(0, 1)], &[1.0]);
        let law = Caloric;
        let old = inventory(&law, &[300.0; 2], &[0.0; 2]);
        let cells = [CellId(0), CellId(1)];
        assert!(
            GraphEnthalpyTransport::from_volume_flows(
                cx,
                &graph,
                &[0.1],
                Density::new(1.0),
                &cells,
                &old,
                vec![],
                3,
            )
            .is_err()
        );
        let binding = GraphEnthalpyTransport::from_volume_flows(
            cx,
            &graph,
            &[0.1],
            Density::new(1.0),
            &cells,
            &old,
            vec![
                MassBoundary::Inlet {
                    id: 10,
                    cell: cells[0],
                    kg_s: 0.1,
                    material_id: MaterialId(1),
                    specific_h_j_kg: 60_000.0,
                },
                MassBoundary::Outlet {
                    id: 11,
                    cell: cells[1],
                    kg_s: 0.1,
                },
            ],
            3,
        )
        .unwrap();
        let (next, audit) = binding
            .advance_advection_implicit(cx, &old, 1.0, 1e-7, 1e-7)
            .unwrap();
        near(
            next.temperature_k(cells[0]).unwrap(),
            300.0 + 60.0 / 11.0,
            1e-8,
        );
        near(
            next.temperature_k(cells[1]).unwrap(),
            300.0 + 60.0 / 121.0,
            1e-8,
        );
        near(audit.stored_energy_change_j, audit.external_heat_j, 1e-7);
    });
}

// A physical two-branch ring with four distinct native pipe inventories.
// Native CellId order deliberately differs from the positive hydraulic route.
const ROUTED_NODES: [CellId; 2] = [CellId(0), CellId(1)];
const ROUTED_ORDER: [CellId; 6] = [
    CellId(0),
    CellId(3),
    CellId(2),
    CellId(1),
    CellId(5),
    CellId(4),
];

fn routed_data() -> Vec<DepositedCell> {
    let geometry = [
        ([0.0, 0.0, 0.0], [0.0, 0.1, 0.0]),
        ([2.0, 0.0, 0.0], [2.0, 0.1, 0.0]),
        ([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]),
        ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
        ([2.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
    ];
    geometry
        .into_iter()
        .enumerate()
        .map(|(i, (start, end))| DepositedCell {
            id: CellId(i as u64),
            segment_id: i as u64,
            body_id: 0,
            material_id: MaterialId(1),
            mass_kg: 1.0,
            specific_h_j_kg: if i == 0 { 60_000.0 } else { 0.0 },
            geometry: CellGeometry {
                start_m: start,
                end_m: end,
                volume_m3: 1.0,
                surface_area_m2: 1.0,
            },
            boundary: Boundary {
                ambient_k: 300.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        })
        .collect()
}

fn routed_state<'a>(law: &'a Caloric, cells: &[DepositedCell]) -> DepositionSimulation<'a> {
    DepositionSimulation::new(
        cells
            .iter()
            .cloned()
            .map(|cell| TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            })
            .collect(),
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: cells.len(),
            max_links: 8,
            max_events: cells.len() + 2,
            max_substeps: 16,
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
}

fn routed_layout(
    cx: &Cx<'_>,
    graph: &LossGraph,
    thermal: &DepositionSimulation<'_>,
) -> DistributedGraphInventory {
    DistributedGraphInventory::new(
        cx,
        graph,
        Density::new(1.0),
        &ROUTED_NODES,
        &[vec![CellId(3), CellId(2)], vec![CellId(5), CellId(4)]],
        thermal,
        16,
    )
    .unwrap()
}

fn routed_momentum(
    cx: &Cx<'_>,
    graph: &LossGraph,
    layout: &DistributedGraphInventory,
    flow: f64,
) -> InertialGraph {
    let sections: Vec<_> = layout
        .inertial_branches()
        .iter()
        .map(|geometry| geometry.expect("every fixture branch has physical pipe inventory"))
        .collect();
    InertialGraph::try_new(
        cx,
        graph,
        layout.reference_density(),
        &sections,
        &[boundary(0, 0.0)],
        &[flow; 2],
        0.0,
        16,
    )
    .unwrap()
}

fn routed_audits<P, T>(report: &InertialThermalReport<P, T>) {
    let limits = controls();
    assert!(report.maximum_cell_residual_j <= limits.cell_tolerance_j);
    assert!(report.thermal.residual_j.abs() <= limits.total_thermal_tolerance_j);
    assert!(report.hydraulics.maximum_impulse_residual_pa_s <= limits.impulse_tolerance_pa_s);
    assert!(report.hydraulics.mechanical_balance_residual_j.abs() <= limits.mechanical_tolerance_j);
    assert!(
        report
            .hydraulics
            .average_flow_residual_m3_s
            .iter()
            .all(|r| r.abs() <= limits.average_flow_tolerance_m3_s)
    );
}

#[test]
fn g1_distributed_pipe_residence_and_reversal_retain_the_original_momentum_dimension() {
    with_cx(|cx| {
        let graph = network(2, &[(0, 1), (1, 0)], &[1.0; 2]);
        let law = Caloric;
        let data = routed_data();
        let old = routed_state(&law, &data);
        let layout = routed_layout(cx, &graph, &old);
        let momentum = routed_momentum(cx, &graph, &layout, 1.0);
        layout.admit_momentum(cx, &momentum).unwrap();
        assert_eq!(layout.fluid_cell_count(), 6);
        assert_eq!(layout.interior_connection_count(), 6);
        assert_eq!(layout.branch_volumes_m3(), &[2.0; 2]);
        assert_eq!(momentum.inertances_pa_s2_m3(), &[2.0; 2]);
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &ROUTED_NODES,
            maximum_connections: 16,
        };
        let (forward_momentum, forward, first) = advance_inertial_graph_distributed(
            cx,
            &problem,
            &layout,
            &old,
            1.0,
            controls(),
            &[0.0; 6],
            &[1.0; 2],
            |_, _, _| Ok::<_, String>(response(&graph, &[1.0; 2])),
            |_, _, _, h| Ok(heat(h.average_volume_flows_m3_s.clone())),
        )
        .unwrap();
        // A closed six-tank BE recurrence with mdot*dt/M=1. Only the first
        // junction starts warm: theta_i=30*(1/2)^i/(1-(1/2)^6).
        let amplitude = 30.0 / (1.0 - 0.5_f64.powi(6));
        for (i, id) in ROUTED_ORDER.into_iter().enumerate() {
            near(
                forward.temperature_k(id).unwrap(),
                300.0 + amplitude * 0.5_f64.powi(i as i32),
                1e-7,
            );
        }
        for q in &first.hydraulics.average_volume_flows_m3_s {
            near(*q, 1.0, 1e-10);
        }
        assert_eq!(first.heat, first.hydraulics.average_volume_flows_m3_s);
        // Thermal refinement introduced no hydraulic node or impulse coordinate.
        assert_eq!(first.hydraulics.momentum.x.len(), 3); // two impulses + one pressure
        assert_eq!(forward_momentum.node_count(), 2);
        assert_eq!(forward_momentum.branch_count(), 2);

        let direct = routed_state(&law, &data[..2]);
        let (_, instant, direct_report) = advance_inertial_graph(
            cx,
            &problem,
            &direct,
            1.0,
            controls(),
            &[0.0; 2],
            &[1.0; 2],
            |_, _, _| Ok::<_, String>(response(&graph, &[1.0; 2])),
            |_, _, _, _| Ok(heat(())),
        )
        .unwrap();
        assert_eq!(
            direct_report.hydraulics.momentum.x.len(),
            first.hydraulics.momentum.x.len()
        );
        near(instant.temperature_k(ROUTED_NODES[1]).unwrap(), 320.0, 1e-7);
        assert!(forward.temperature_k(ROUTED_NODES[1]).unwrap() < 304.0);

        let reversal_problem = InertialThermalProblem {
            momentum: &forward_momentum,
            node_cells: &ROUTED_NODES,
            maximum_connections: 16,
        };
        let (back_momentum, back, reversed) = advance_inertial_graph_distributed(
            cx,
            &reversal_problem,
            &layout,
            &forward,
            2.0,
            controls(),
            &[0.0; 6],
            &[-1.0; 2],
            |_, _, _| Ok::<_, String>(response(&graph, &[-9.0; 2])),
            |_, _, _, h| Ok(heat(h.average_volume_flows_m3_s.clone())),
        )
        .unwrap();
        // Native midpoint: (2I/dt+R)*Qavg=2I*Qold/dt+drive -> Qavg=-1,
        // Qnew=-3. This checks the documented average-flow discretization;
        // it does not claim to resolve the within-step zero-crossing time.
        for edge in 0..2 {
            near(
                reversed.hydraulics.average_volume_flows_m3_s[edge],
                -1.0,
                1e-10,
            );
            near(reversed.hydraulics.final_flows_m3_s[edge], -3.0, 1e-10);
        }
        // Reverse the SAME stored spatial profile. The independent periodic
        // recurrence now receives its donor from the next positive-route cell.
        for (i, id) in ROUTED_ORDER.into_iter().enumerate() {
            let weighted: f64 = (0..6)
                .map(|k| {
                    let donor = ROUTED_ORDER[(i + k) % 6];
                    (forward.temperature_k(donor).unwrap() - 300.0) * 0.5_f64.powi(k as i32 + 1)
                })
                .sum();
            let expected = 300.0 + weighted / (1.0 - 0.5_f64.powi(6));
            near(back.temperature_k(id).unwrap(), expected, 1e-7);
        }
        assert!(back.temperature_k(CellId(4)).unwrap() > forward.temperature_k(CellId(4)).unwrap());
        assert_eq!(reversed.hydraulics.momentum.x.len(), 3);
        assert_eq!(
            back_momentum.inertances_pa_s2_m3(),
            momentum.inertances_pa_s2_m3()
        );
        assert_eq!(back_momentum.accepted_time_s(), 2.0);
        assert_eq!(back.time_s(), 2.0);
        for (&id, original) in old.cells() {
            let retained = &back.cells()[&id];
            assert_eq!(retained.mass_kg, original.mass_kg);
            assert_eq!(retained.geometry, original.geometry);
            assert_eq!(retained.material_id, original.material_id);
        }
        for report in [&first, &reversed] {
            routed_audits(report);
            near(report.thermal.external_heat_j, 0.0, 1e-12);
            near(report.thermal.stored_energy_change_j, 0.0, 1e-6);
        }
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(momentum.accepted_time_s(), 0.0);
        assert_eq!(momentum.volume_flows_m3_s(), vec![1.0; 2]);
    });
}

#[test]
fn g1_distributed_pipe_heat_changes_drag_at_the_accepted_native_endpoint() {
    with_cx(|cx| {
        let graph = network(2, &[(0, 1), (1, 0)], &[1.0; 2]);
        let law = Caloric;
        let mut data = routed_data();
        data[0].specific_h_j_kg = 0.0;
        data[3].boundary.source_w = 1000.0;
        let old = routed_state(&law, &data);
        let layout = routed_layout(cx, &graph, &old);
        let momentum = routed_momentum(cx, &graph, &layout, 0.0);
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &ROUTED_NODES,
            maximum_connections: 16,
        };
        let (_, next, report) = advance_inertial_graph_distributed(
            cx,
            &problem,
            &layout,
            &old,
            1.0,
            controls(),
            &[0.0; 6],
            &[0.0; 2],
            |_, candidate, _| {
                let t = candidate
                    .temperature_k(CellId(3))
                    .map_err(|e| e.to_string())?;
                let r = 1.0 + 0.1 * (t - 300.0);
                Ok::<_, String>(ResolvedInertialGraph {
                    graph: network(2, &[(0, 1), (1, 0)], &[r, 1.0]),
                    pressure_boundaries: vec![boundary(0, 0.0)],
                    driving_pressure_pa: vec![1.0; 2],
                    payload: r,
                })
            },
            |_, _, &r, _| Ok(heat(r)),
        )
        .unwrap();
        // Eliminate pressure around the two-branch ring. Native inertance=2
        // each: Q*(8 + R_hot + 1)=2. The six-cell periodic heat recurrence
        // gives theta_hot=1/[(1+Q)*(1-(Q/(1+Q))^6)].
        let hot = |q: f64| 1.0 / ((1.0 + q) * (1.0 - (q / (1.0 + q)).powi(6)));
        let (mut low, mut high) = (0.0, 0.2);
        for _ in 0..80 {
            let q = 0.5 * (low + high);
            if q * (10.0 + 0.1 * hot(q)) > 2.0 {
                high = q;
            } else {
                low = q;
            }
        }
        let q = 0.5 * (low + high);
        near(next.temperature_k(CellId(3)).unwrap(), 300.0 + hot(q), 1e-7);
        for rate in &report.hydraulics.average_volume_flows_m3_s {
            near(*rate, q, 1e-9);
        }
        for rate in &report.hydraulics.final_flows_m3_s {
            near(*rate, 2.0 * q, 1e-9);
        }
        near(report.sources, 1.0 + 0.1 * hot(q), 1e-9);
        near(
            report.sources,
            1.0 + 0.1 * (next.temperature_k(CellId(3)).unwrap() - 300.0),
            1e-12,
        );
        assert_eq!(report.heat, report.sources);
        assert!(q < 0.2);
        near(report.thermal.external_heat_j, 1000.0, 1e-8);
        near(report.thermal.stored_energy_change_j, 1000.0, 1e-6);
        routed_audits(&report);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(momentum.accepted_time_s(), 0.0);
    });
}

#[test]
fn g0_distributed_momentum_and_storage_mismatches_publish_neither_owner() {
    with_cx(|cx| {
        let graph = network(2, &[(0, 1), (1, 0)], &[1.0; 2]);
        let law = Caloric;
        let data = routed_data();
        let old = routed_state(&law, &data);
        let layout = routed_layout(cx, &graph, &old);
        let momentum = routed_momentum(cx, &graph, &layout, 1.0);
        let wrong = owner(cx, &graph, &[3.0; 2], &[1.0; 2]);
        let calls = Cell::new(0usize);
        let wrong_problem = InertialThermalProblem {
            momentum: &wrong,
            node_cells: &ROUTED_NODES,
            maximum_connections: 16,
        };
        let result = advance_inertial_graph_distributed(
            cx,
            &wrong_problem,
            &layout,
            &old,
            1.0,
            controls(),
            &[0.0; 6],
            &[1.0; 2],
            |_, _, _| {
                calls.set(calls.get() + 1);
                Ok::<_, String>(response(&graph, &[1.0; 2]))
            },
            |_, _, _, _| Ok(heat(())),
        );
        assert!(matches!(result, Err(InertialThermalError::Transport(_))));
        assert_eq!(calls.get(), 0); // wrong physical inertia refuses before source work

        let mut changed = data.clone();
        changed[3].geometry.end_m[0] += 0.25;
        let changed = routed_state(&law, &changed);
        let problem = InertialThermalProblem {
            momentum: &momentum,
            node_cells: &ROUTED_NODES,
            maximum_connections: 16,
        };
        let result = advance_inertial_graph_distributed(
            cx,
            &problem,
            &layout,
            &changed,
            1.0,
            controls(),
            &[0.0; 6],
            &[1.0; 2],
            |_, _, _| {
                calls.set(calls.get() + 1);
                Ok::<_, String>(response(&graph, &[1.0; 2]))
            },
            |_, _, _, _| Ok(heat(())),
        );
        assert!(matches!(result, Err(InertialThermalError::Transport(_))));
        assert!(calls.get() > 0); // original candidate storage admission rejected the change
        for state in [&old, &changed] {
            assert_eq!(state.time_s(), 0.0);
        }
        for owner in [&momentum, &wrong] {
            assert_eq!(owner.accepted_time_s(), 0.0);
            assert_eq!(owner.volume_flows_m3_s(), vec![1.0; 2]);
        }
        for (&id, cell) in old.cells() {
            assert_eq!(cell.specific_h_j_kg, data[id.0 as usize].specific_h_j_kg);
        }
    });
}
