//! Actual material/hydraulic/thermal owners against independent small-system roots.

use super::*;
use crate::graph::porous::{LinearLossElement, LinearLossResistance, PassiveBranch};
use crate::graph::tests::{boundary, config as graph_controls, edge, with_cx};
use crate::{SourceProvenance, ToleranceBasis};
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
