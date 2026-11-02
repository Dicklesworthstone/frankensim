//! Independent finite-network roots exercise the actual momentum and heat owners.
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
