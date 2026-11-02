//! Actual material/hydraulic/thermal owners against independent small-system roots.

use super::*;
use crate::graph::gravity::GraphGravitySection;
use crate::graph::porous::{LinearLossElement, LinearLossResistance, PassiveBranch};
use crate::graph::tests::{boundary, config as graph_controls, edge, with_cx};
use crate::{SourceProvenance, ToleranceBasis};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::gas::{GasSpec, GasState};
use fs_qty::{Length, VolumetricFlowRate};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, MaterialId, StepConfig,
    ThermalLaw, TimedEvent,
};
use std::collections::BTreeMap;

struct GasCaloric {
    cp: f64,
}
impl ThermalLaw for GasCaloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h / self.cp;
        if t.is_finite() && (200.0..=500.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp
    }
}
fn gas(t: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), t, 101325.0).unwrap()
}
fn inventory<'a>(
    law: &'a GasCaloric,
    rho: f64,
    t: [f64; 2],
    powers: [f64; 2],
) -> DepositionSimulation<'a> {
    let schedule = (0..2)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: 0,
                material_id: MaterialId(1),
                mass_kg: 1.0,
                specific_h_j_kg: law.cp * (t[i] - 300.0),
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0 / rho,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 400.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: powers[i],
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 8,
            max_events: 4,
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
fn controls() -> QuasiSteadyGraphControls {
    QuasiSteadyGraphControls {
        solver: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        maximum_cells: 2,
        maximum_newton_iterations: 30,
        maximum_endpoint_evaluations: 512,
        energy_scale_j: 1000.0,
        finite_difference_step: 1e-5,
        cell_tolerance_j: 1e-6,
        total_tolerance_j: 1e-6,
        maximum_neglected_power_w: 100.0,
        maximum_neglected_energy_j: 1000.0,
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
fn g0_quasi_steady_work_counts_opposing_local_heads_before_they_cancel() {
    with_cx(|cx| {
        let material = |rho| FluidTransportState::try_new(rho, 1000.0, 0.1, 0.001).unwrap();
        let points = [[0.0; 3], [0.0, 0.0, 2.0], [1.0, 0.0, 2.0], [1.0, 0.0, 0.0]];
        let sections = (0..3)
            .map(|i| GraphGravitySection {
                start_m: points[i],
                end_m: points[i + 1],
                transport: material(0.75),
            })
            .collect();
        let graph = LossGraph::new(2, vec![edge("bent", 0, 1, 1.0)])
            .unwrap()
            .with_distributed_boussinesq_gravity(
                &material(1.0),
                &[sections],
                &[points[0], points[3]],
                [0.0, 0.0, -8.0],
                0.3,
                3,
                cx,
            )
            .unwrap();
        let flow = graph
            .solve(&[boundary(0, 1.0), boundary(1, 0.0)], graph_controls(), cx)
            .unwrap();
        assert_eq!(flow.branches[0].gravity_pressure_head.value(), 0.0);
        assert_eq!(flow.branches[0].flow.value(), 1.0);
        // 1 W drag plus |+4 Pa|*1 m3/s + |-4 Pa|*1 m3/s.
        assert_eq!(neglected_power::<()>(cx, &flow).unwrap(), 9.0);
    });
}

const PIPE_CELL_MASS_KG: f64 = 0.001;

fn pipe_reference_resistance() -> f64 {
    let reference = gas(300.0);
    // Two circular, one-metre cells: A=V/L and Hagen-Poiseuille R=8*pi*mu*L/A^2.
    let area = PIPE_CELL_MASS_KG / reference.density;
    16.0 * core::f64::consts::PI * reference.dynamic_viscosity / (area * area)
}

fn pipe_inventory<'a>(
    law: &'a dyn ThermalLaw,
    rho: f64,
    middle_x: f64,
) -> DepositionSimulation<'a> {
    // The hydraulic graph has only nodes 0 and 1. Pipe order is deliberately
    // 20 then 10, so the authored path cannot be replaced by CellId sorting.
    let geometry = [
        (0, [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
        (1, [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]),
        (20, [0.0, 0.0, 0.0], [middle_x, 0.0, 0.0]),
        (10, [middle_x, 0.0, 0.0], [2.0, 0.0, 0.0]),
    ];
    let schedule = geometry
        .into_iter()
        .map(|(id, start_m, end_m)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: 0,
                material_id: MaterialId(1),
                mass_kg: PIPE_CELL_MASS_KG,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m,
                    end_m,
                    volume_m3: PIPE_CELL_MASS_KG / rho,
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
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 4,
            max_links: 16,
            max_events: 4,
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

fn distributed_response(
    cx: &Cx<'_>,
    candidate: &DepositionSimulation<'_>,
    sign: f64,
    inlet_excess_k: f64,
) -> Result<QuasiSteadyGraphResponse<GasState>, String> {
    let mean_pipe_temperature = 0.5
        * (candidate
            .temperature_k(CellId(20))
            .map_err(|e| e.to_string())?
            + candidate
                .temperature_k(CellId(10))
                .map_err(|e| e.to_string())?);
    let reference = gas(300.0);
    let source = gas(mean_pipe_temperature);
    let graph = viscous_graph(
        pipe_reference_resistance() * source.dynamic_viscosity / reference.dynamic_viscosity,
        false,
    );
    let hydraulics = graph
        .solve(
            &[
                boundary(
                    0,
                    sign * pipe_reference_resistance() * PIPE_CELL_MASS_KG / reference.density,
                ),
                boundary(1, 0.0),
            ],
            graph_controls(),
            cx,
        )
        .map_err(|e| e.to_string())?;
    let mass = reference.density * hydraulics.branches[0].flow.value().abs();
    let area = PIPE_CELL_MASS_KG / reference.density;
    let diameter = (4.0 * area / core::f64::consts::PI).sqrt();
    let reynolds = mass * diameter / (area * source.dynamic_viscosity);
    if !(0.0..=2300.0).contains(&reynolds) {
        return Err("distributed capillary fixture left its laminar source support".to_string());
    }
    let inlet = if sign > 0.0 { CellId(0) } else { CellId(1) };
    let outlet = if sign > 0.0 { CellId(1) } else { CellId(0) };
    Ok(QuasiSteadyGraphResponse {
        hydraulics,
        mass_boundaries: vec![
            MassBoundary::Inlet {
                id: 100,
                cell: inlet,
                kg_s: mass,
                material_id: MaterialId(1),
                specific_h_j_kg: reference.specific_heat_cp * inlet_excess_k,
            },
            MassBoundary::Outlet {
                id: 101,
                cell: outlet,
                kg_s: mass,
            },
        ],
        heat: heat(source),
    })
}

#[test]
fn g1_distributed_quasi_steady_residence_and_reversal_keep_the_original_hydraulics() {
    with_cx(|cx| {
        let reference = gas(300.0);
        let law = GasCaloric {
            cp: reference.specific_heat_cp,
        };
        let graph = viscous_graph(pipe_reference_resistance(), false);
        let nodes = [CellId(0), CellId(1)];
        let old = pipe_inventory(&law, reference.density, 1.0);
        let inventory = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(reference.density),
            &nodes,
            &[vec![CellId(20), CellId(10)]],
            &old,
            8,
        )
        .unwrap();
        let coupling = problem(&graph, &nodes, reference.density);
        let mut limits = controls();
        limits.maximum_cells = 4;
        limits.energy_scale_j = 1.0;
        let ids = [CellId(0), CellId(1), CellId(20), CellId(10)];
        let mut before = [0.0; 4];
        let mut current = old.clone();
        for (step_index, (sign, inlet)) in [(1.0, 60.0), (-1.0, 0.0)].into_iter().enumerate() {
            // With a=mdot*dt/M, native BE is the ordered recurrence
            // e_new=(e_old+a*e_upstream)/(1+a). Solve only the capillary relation
            // a*mu(T_pipe_mean)/mu(300K)=1 independently of production Newton.
            let order = if sign > 0.0 {
                [0, 2, 3, 1]
            } else {
                [1, 3, 2, 0]
            };
            let field = |turnover: f64| {
                let mut next = [0.0; 4];
                let mut upstream = inlet;
                for index in order {
                    next[index] = (before[index] + turnover * upstream) / (1.0 + turnover);
                    upstream = next[index];
                }
                next
            };
            let (mut lo, mut hi) = (0.0, 1.0);
            for _ in 0..80 {
                let turnover = 0.5 * (lo + hi);
                let temperatures = field(turnover);
                let resistance = gas(300.0 + 0.5 * (temperatures[2] + temperatures[3]))
                    .dynamic_viscosity
                    / reference.dynamic_viscosity;
                if turnover * resistance > 1.0 {
                    hi = turnover;
                } else {
                    lo = turnover;
                }
            }
            let turnover = 0.5 * (lo + hi);
            let expected = field(turnover);
            let (next, report) = advance_quasi_steady_graph_distributed(
                cx,
                &coupling,
                &inventory,
                &current,
                (step_index + 1) as f64,
                limits,
                &[0.0; 4],
                |cx, candidate| distributed_response(cx, candidate, sign, inlet),
            )
            .unwrap();
            assert_eq!(report.hydraulics.pressures.len(), 2);
            assert_eq!(report.hydraulics.branches.len(), 1);
            assert_eq!(next.cells().len(), 4);
            let q = report.hydraulics.branches[0].flow.value();
            near(
                q,
                sign * PIPE_CELL_MASS_KG * turnover / reference.density,
                2e-12,
            );
            for (i, id) in ids.into_iter().enumerate() {
                near(next.temperature_k(id).unwrap(), 300.0 + expected[i], 2e-7);
            }
            near(
                report.heat.dynamic_viscosity,
                gas(300.0 + 0.5 * (expected[2] + expected[3])).dynamic_viscosity,
                1e-12,
            );
            let residence = inventory.nominal_residence_times_s(cx, &[q]).unwrap()[0].unwrap();
            near(residence, 2.0 / turnover, 2e-8);
            let expected_energy = PIPE_CELL_MASS_KG
                * law.cp
                * (expected.iter().sum::<f64>() - before.iter().sum::<f64>());
            near(report.thermal.stored_energy_change_j, expected_energy, 1e-6);
            near(report.thermal.external_heat_j, expected_energy, 1e-6);
            if sign < 0.0 {
                assert!(report.thermal.external_heat_j < 0.0);
            }
            before = expected;
            current = next;
        }
        assert_eq!(old.time_s(), 0.0);
        assert!(old.cells().values().all(|cell| cell.specific_h_j_kg == 0.0));
        assert_eq!(current.time_s(), 2.0);
    });
}

#[test]
fn g0_distributed_quasi_steady_admission_preserves_inventory_and_source_identity() {
    struct TaggedCaloric {
        law: GasCaloric,
        identity: u8,
    }
    impl ThermalLaw for TaggedCaloric {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            self.law.temperature_k(h)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            self.law.maximum_temperature_slope()
        }
        fn caloric_identity(&self) -> Option<[u8; 32]> {
            Some([self.identity; 32])
        }
    }
    with_cx(|cx| {
        let reference = gas(300.0);
        let law = TaggedCaloric {
            law: GasCaloric {
                cp: reference.specific_heat_cp,
            },
            identity: 1,
        };
        let other_law = TaggedCaloric {
            law: GasCaloric {
                cp: reference.specific_heat_cp,
            },
            identity: 2,
        };
        let graph = viscous_graph(pipe_reference_resistance(), false);
        let reversed_graph = viscous_graph(pipe_reference_resistance(), true);
        let nodes = [CellId(0), CellId(1)];
        let swapped = [CellId(1), CellId(0)];
        let old = pipe_inventory(&law, reference.density, 1.0);
        let changed_geometry = pipe_inventory(&law, reference.density, 1.25);
        let changed_source = pipe_inventory(&other_law, reference.density, 1.0);
        let inventory = DistributedGraphInventory::new(
            cx,
            &graph,
            Density::new(reference.density),
            &nodes,
            &[vec![CellId(20), CellId(10)]],
            &old,
            8,
        )
        .unwrap();
        let good = problem(&graph, &nodes, reference.density);
        let wrong_nodes = problem(&graph, &swapped, reference.density);
        let wrong_density = problem(&graph, &nodes, 2.0 * reference.density);
        let wrong_topology = problem(&reversed_graph, &nodes, reference.density);
        let mut insufficient_connections = problem(&graph, &nodes, reference.density);
        insufficient_connections.maximum_connections = 2;
        let mut limits = controls();
        limits.maximum_cells = 4;
        for (problem, state) in [
            (&wrong_nodes, &old),
            (&wrong_density, &old),
            (&wrong_topology, &old),
            (&insufficient_connections, &old),
            (&good, &changed_geometry),
            (&good, &changed_source),
        ] {
            let refused = advance_quasi_steady_graph_distributed(
                cx,
                problem,
                &inventory,
                state,
                1.0,
                limits,
                &[0.0; 4],
                |_, _| -> Result<QuasiSteadyGraphResponse<()>, &str> {
                    panic!("inventory mismatch must refuse before source callback work")
                },
            );
            assert!(matches!(
                refused,
                Err(QuasiSteadyGraphError::Invalid(_)) | Err(QuasiSteadyGraphError::Transport(_))
            ));
        }
        // A solved response also retains its actual gravity chart. The base
        // graph has no gravity, so checking only that template would miss this.
        let wrong_chart = graph
            .clone()
            .with_gravity(
                &reference,
                &[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
                [0.0; 3],
                cx,
            )
            .unwrap();
        let refused = advance_quasi_steady_graph_distributed(
            cx,
            &good,
            &inventory,
            &old,
            1.0,
            limits,
            &[0.0; 4],
            |cx, candidate| {
                let mut response = distributed_response(cx, candidate, 1.0, 60.0)?;
                response.hydraulics = wrong_chart
                    .solve(
                        &[
                            boundary(
                                0,
                                pipe_reference_resistance() * PIPE_CELL_MASS_KG / reference.density,
                            ),
                            boundary(1, 0.0),
                        ],
                        graph_controls(),
                        cx,
                    )
                    .map_err(|e| e.to_string())?;
                Ok::<_, String>(response)
            },
        );
        assert!(matches!(
            refused,
            Err(QuasiSteadyGraphError::Transport(
                TransientTransportError::Invalid(
                    "hydraulic gravity and distributed reference inventory disagree"
                )
            ))
        ));
        assert_eq!(old.time_s(), 0.0);
    });
}

fn problem<'a>(graph: &'a LossGraph, cells: &'a [CellId], rho: f64) -> QuasiSteadyGraphProblem<'a> {
    QuasiSteadyGraphProblem {
        topology: graph,
        node_cells: cells,
        reference_density: Density::new(rho),
        maximum_connections: 8,
    }
}
fn loop_response(
    cx: &Cx<'_>,
    graph: &LossGraph,
    candidate: &DepositionSimulation<'_>,
) -> Result<QuasiSteadyGraphResponse<[f64; 2]>, String> {
    let temperatures = [
        candidate
            .temperature_k(CellId(0))
            .map_err(|e| e.to_string())?,
        candidate
            .temperature_k(CellId(1))
            .map_err(|e| e.to_string())?,
    ];
    let materials = [gas(temperatures[0]), gas(temperatures[1])];
    let mut settings = graph_controls();
    settings.absolute_flow_tolerance = VolumetricFlowRate::new(1e-17);
    settings.relative_flow_tolerance = 1e-15;
    let hydraulics = graph
        .clone()
        .with_boussinesq_gravity(
            &gas(315.0),
            &materials,
            &[[0.0; 3], [0.0, 0.0, 1.0]],
            [0.0, 0.0, -9.81],
            0.25,
            cx,
        )
        .map_err(|e| e.to_string())?
        .solve(&[boundary(0, 0.0)], settings, cx)
        .map_err(|e| e.to_string())?;
    Ok(QuasiSteadyGraphResponse {
        hydraulics,
        mass_boundaries: vec![],
        heat: heat(temperatures),
    })
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn g1_endpoint_material_buoyancy_and_enthalpy_match_the_implicit_thermosiphon_root() {
    with_cx(|cx| {
        let reference = gas(315.0);
        let law = GasCaloric {
            cp: reference.specific_heat_cp,
        };
        let graph = LossGraph::new(
            2,
            vec![edge("rise", 0, 1, 10_000.0), edge("return", 1, 0, 10_000.0)],
        )
        .unwrap();
        let cells = [CellId(0), CellId(1)];
        let coupling = problem(&graph, &cells, reference.density);
        for sign in [-1.0, 1.0] {
            let old = inventory(
                &law,
                reference.density,
                [315.0 + 15.0 * sign, 315.0 - 15.0 * sign],
                [1000.0 * sign, -1000.0 * sign],
            );
            let initial_q = loop_response(cx, &graph, &old).unwrap().hydraulics.branches[0]
                .flow
                .value();
            let (next, report) = advance_quasi_steady_graph(
                cx,
                &coupling,
                &old,
                10.0,
                controls(),
                &[0.0; 2],
                |cx, state| loop_response(cx, &graph, state),
            )
            .unwrap();
            // Equal inventories: dT = (dT_old + 2 dt P/Cp) / (1 + 2 dt rho0 |Q|).
            // Ideal-gas buoyancy supplies 2 R Q^2 = g (rho_cold-rho_hot).
            // Solve this independent scalar equation without production iteration.
            let numerator = 30.0 + 20_000.0 / law.cp;
            let flow_for = |difference: f64| {
                (9.81
                    * reference.density
                    * 315.0
                    * (1.0 / (315.0 - 0.5 * difference) - 1.0 / (315.0 + 0.5 * difference))
                    / 20_000.0)
                    .sqrt()
            };
            let (mut lo, mut hi) = (0.0, numerator);
            for _ in 0..80 {
                let d = 0.5 * (lo + hi);
                if d * (1.0 + 20.0 * reference.density * flow_for(d)) > numerator {
                    hi = d;
                } else {
                    lo = d;
                }
            }
            let difference = 0.5 * (lo + hi);
            let q = report.hydraulics.branches[0].flow.value();
            near(q, sign * flow_for(difference), 2e-10);
            near(
                next.temperature_k(CellId(0)).unwrap(),
                315.0 + sign * 0.5 * difference,
                2e-7,
            );
            near(
                next.temperature_k(CellId(1)).unwrap(),
                315.0 - sign * 0.5 * difference,
                2e-7,
            );
            assert!(q.abs() > initial_q.abs());
            near(report.thermal.stored_energy_change_j, 0.0, 1e-6);
            near(report.thermal.external_heat_j, 0.0, 1e-6);
            near(
                report.heat[0],
                next.temperature_k(CellId(0)).unwrap(),
                1e-12,
            );
            assert_eq!(old.time_s(), 0.0);
            assert_eq!(next.time_s(), 10.0);
        }
    });
}

fn viscous_graph(resistance: f64, reversed: bool) -> LossGraph {
    let loss = LinearLossElement::new(
        "duct",
        LinearLossResistance::new(resistance),
        0.0,
        SourceProvenance::new(
            "material viscosity and authored capillary geometry",
            "quasi-steady-test-v1",
        ),
        ToleranceBasis::EngineeringAllowance,
    )
    .unwrap();
    LossGraph::new_mixed(
        2,
        vec![PassiveBranch::Linear {
            from: usize::from(reversed),
            to: usize::from(!reversed),
            loss,
        }],
    )
    .unwrap()
}
fn pressure_response(
    cx: &Cx<'_>,
    candidate: &DepositionSimulation<'_>,
    wall_exchange: bool,
    reversed: bool,
) -> Result<QuasiSteadyGraphResponse<[f64; 2]>, String> {
    let temperatures = [
        candidate
            .temperature_k(CellId(0))
            .map_err(|e| e.to_string())?,
        candidate
            .temperature_k(CellId(1))
            .map_err(|e| e.to_string())?,
    ];
    let source = gas(0.5 * (temperatures[0] + temperatures[1]));
    let graph = viscous_graph(
        source.dynamic_viscosity / gas(300.0).dynamic_viscosity,
        reversed,
    );
    let hydraulics = graph
        .solve(&[boundary(0, 1.0), boundary(1, 0.0)], graph_controls(), cx)
        .map_err(|e| e.to_string())?;
    let mass = gas(300.0).density * hydraulics.branches[0].flow.value().abs();
    let mut heat = heat([0.0; 2]);
    if wall_exchange {
        use crate::graph::thermal::medium::convection::{
            ForcedConvectionSettings, evaluate_forced_from_material,
        };
        for (index, &temperature) in temperatures.iter().enumerate() {
            let material = gas(temperature);
            let exchange = evaluate_forced_from_material(
                fs_convection::CorrelationId::CircularDuctHausen,
                &material,
                hydraulics.branches[0].flow.value().abs(),
                Length::new(0.01),
                ForcedConvectionSettings {
                    length_over_hydraulic_diameter: Some(100.0),
                    ..ForcedConvectionSettings::default()
                },
            )
            .map_err(|e| e.to_string())?;
            let coefficient = exchange.coefficient().value.value();
            heat.boundaries
                .ambient_w_k
                .insert(CellId(index as u64), coefficient);
            heat.bulk.payload[index] = coefficient;
        }
    }
    Ok(QuasiSteadyGraphResponse {
        hydraulics,
        mass_boundaries: vec![
            MassBoundary::Inlet {
                id: 20,
                cell: CellId(0),
                kg_s: mass,
                material_id: MaterialId(1),
                specific_h_j_kg: gas(300.0).specific_heat_cp * 60.0,
            },
            MassBoundary::Outlet {
                id: 21,
                cell: CellId(1),
                kg_s: mass,
            },
        ],
        heat,
    })
}

#[test]
fn g1_source_viscosity_closes_open_transport_with_reversed_edge_declarations() {
    with_cx(|cx| {
        let reference = gas(300.0);
        let law = GasCaloric {
            cp: reference.specific_heat_cp,
        };
        let old = inventory(&law, reference.density, [300.0; 2], [0.0; 2]);
        let cells = [CellId(0), CellId(1)];
        for reversed in [false, true] {
            let graph = viscous_graph(1.0, reversed);
            let coupling = problem(&graph, &cells, reference.density);
            let (next, report) = advance_quasi_steady_graph(
                cx,
                &coupling,
                &old,
                1.0,
                controls(),
                &[0.0; 2],
                |cx, state| pressure_response(cx, state, false, reversed),
            )
            .unwrap();
            // Two equal-volume CSTR backward-Euler equations, plus the original
            // source's Sutherland viscosity ratio at their mean temperature.
            let temperature_for = |q: f64| {
                let a = reference.density * q / (1.0 + reference.density * q);
                [300.0 + 60.0 * a, 300.0 + 60.0 * a * a]
            };
            let ratio = |t: f64| {
                let s = GasSpec::dry_air_ussa1976().sutherland_s;
                (t / 300.0).powf(1.5) * (300.0 + s) / (t + s)
            };
            let (mut lo, mut hi) = (0.0, 1.0);
            for _ in 0..80 {
                let q = 0.5 * (lo + hi);
                let t = temperature_for(q);
                if q * ratio(0.5 * (t[0] + t[1])) > 1.0 {
                    hi = q;
                } else {
                    lo = q;
                }
            }
            let q = 0.5 * (lo + hi);
            let expected = temperature_for(q);
            near(report.hydraulics.branches[0].flow.value().abs(), q, 1e-9);
            assert!(q < 1.0);
            for (i, t) in expected.into_iter().enumerate() {
                near(next.temperature_k(CellId(i as u64)).unwrap(), t, 2e-7);
            }
            assert!(report.thermal.external_heat_j > 0.0);
            near(
                report.thermal.stored_energy_change_j,
                report.thermal.external_heat_j,
                1e-6,
            );
        }
    });
}

#[test]
fn g3_actual_flow_and_material_update_native_wall_coefficients_in_the_same_step() {
    with_cx(|cx| {
        let reference = gas(300.0);
        let law = GasCaloric {
            cp: reference.specific_heat_cp,
        };
        let old = inventory(&law, reference.density, [300.0; 2], [0.0; 2]);
        let graph = viscous_graph(1.0, false);
        let cells = [CellId(0), CellId(1)];
        let coupling = problem(&graph, &cells, reference.density);
        let (plain, plain_report) = advance_quasi_steady_graph(
            cx,
            &coupling,
            &old,
            1.0,
            controls(),
            &[0.0; 2],
            |cx, state| pressure_response(cx, state, false, false),
        )
        .unwrap();
        let (heated, report) = advance_quasi_steady_graph(
            cx,
            &coupling,
            &old,
            1.0,
            controls(),
            &[0.0; 2],
            |cx, state| pressure_response(cx, state, true, false),
        )
        .unwrap();
        for id in [CellId(0), CellId(1)] {
            assert!(heated.temperature_k(id).unwrap() > plain.temperature_k(id).unwrap());
        }
        assert!(
            report.hydraulics.branches[0].flow.value()
                < plain_report.hydraulics.branches[0].flow.value()
        );
        let independently_resolved = pressure_response(cx, &heated, true, false).unwrap();
        for (actual, expected) in report
            .heat
            .iter()
            .zip(independently_resolved.heat.bulk.payload)
        {
            near(*actual, expected, 1e-10);
        }
        assert!(report.maximum_cell_residual_j <= controls().cell_tolerance_j);
        assert!(report.thermal.residual_j.abs() <= controls().total_tolerance_j);
    });
}

#[test]
fn g0_refusals_leave_inventory_and_time_untouched() {
    with_cx(|cx| {
        let reference = gas(300.0);
        let law = GasCaloric {
            cp: reference.specific_heat_cp,
        };
        let old = inventory(&law, reference.density, [300.0; 2], [0.0; 2]);
        let graph = viscous_graph(1.0, false);
        let cells = [CellId(0), CellId(1)];
        let coupling = problem(&graph, &cells, reference.density);
        let mut limited = controls();
        limited.maximum_endpoint_evaluations = 1;
        assert!(matches!(
            advance_quasi_steady_graph(
                cx,
                &coupling,
                &old,
                1.0,
                limited,
                &[0.0; 2],
                |cx, state| pressure_response(cx, state, false, false)
            ),
            Err(QuasiSteadyGraphError::WorkBudget { .. })
        ));
        let mut limited = controls();
        limited.maximum_neglected_power_w = 0.0;
        assert!(matches!(
            advance_quasi_steady_graph(
                cx,
                &coupling,
                &old,
                1.0,
                limited,
                &[0.0; 2],
                |cx, state| pressure_response(cx, state, false, false)
            ),
            Err(QuasiSteadyGraphError::NeglectedWork { .. })
        ));
        assert!(matches!(
            advance_quasi_steady_graph(
                cx,
                &coupling,
                &old,
                1.0,
                controls(),
                &[0.0; 2],
                |_, _| Err::<QuasiSteadyGraphResponse<()>, _>("source domain")
            ),
            Err(QuasiSteadyGraphError::Resolve("source domain"))
        ));
        assert!(matches!(
            advance_quasi_steady_graph(
                cx,
                &coupling,
                &old,
                1.0,
                controls(),
                &[0.0; 2],
                |cx, state| pressure_response(cx, state, false, true)
            ),
            Err(QuasiSteadyGraphError::Invalid(_))
        ));
        assert_eq!(old.time_s(), 0.0);
        for cell in old.cells().values() {
            assert_eq!(cell.specific_h_j_kg, 0.0);
        }
    });
}

#[test]
fn g4_cancellation_is_before_any_endpoint_publication() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let gate = CancelGate::new_clock_free();
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
        let reference = gas(300.0);
        let law = GasCaloric {
            cp: reference.specific_heat_cp,
        };
        let old = inventory(&law, reference.density, [300.0; 2], [0.0; 2]);
        let graph = viscous_graph(1.0, false);
        let cells = [CellId(0), CellId(1)];
        let coupling = problem(&graph, &cells, reference.density);
        gate.request();
        assert!(matches!(
            advance_quasi_steady_graph(
                &cx,
                &coupling,
                &old,
                1.0,
                controls(),
                &[0.0; 2],
                |_, _| -> Result<QuasiSteadyGraphResponse<()>, String> {
                    panic!("cancelled endpoint must not resolve")
                }
            ),
            Err(QuasiSteadyGraphError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(old.time_s(), 0.0);
    });
}
