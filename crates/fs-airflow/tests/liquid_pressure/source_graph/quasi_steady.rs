//! Original liquid sources at solved quasi-steady thermal endpoints.

use super::*;
use fs_airflow::graph::GraphSolution;
use fs_airflow::graph::thermal::transport::transient::coupled::{
    QuasiSteadyGraphControls, QuasiSteadyGraphProblem, QuasiSteadyGraphResponse,
    advance_quasi_steady_graph_distributed,
};
use fs_solver::NewtonKrylovConfig;

fn hydraulic_controls() -> GraphSolveConfig {
    GraphSolveConfig {
        max_sweeps: 256,
        max_node_iterations: 96,
        absolute_flow_tolerance: VolumetricFlowRate::new(1e-14),
        relative_flow_tolerance: 1e-12,
    }
}

fn fixed_solution(cx: &Cx<'_>, graph: &LossGraph, pressures: [f64; 2]) -> GraphSolution {
    graph
        .solve(
            &pressures
                .into_iter()
                .enumerate()
                .map(|(node, pressure)| FixedPressure {
                    node,
                    pressure: Pressure::new(pressure),
                })
                .collect::<Vec<_>>(),
            hydraulic_controls(),
            cx,
        )
        .unwrap()
}

#[test]
fn bent_source_pipe_advances_native_quasi_steady_buoyancy_and_reciprocal_wall_heat() {
    context(|cx| {
        let card = buoyant_card();
        let (law, caloric) = source_pair(&card);
        // The second pipe cell starts warm and drives the closed loop. Wall
        // heat enters the first cell through its existing reciprocal contact.
        let initial = thermal(
            &caloric,
            [300.0, 300.0, 300.0, 302.0, 300.0, 300.0],
            Some(400.0),
            true,
        );
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let mut laws = laws(2300.0);
        let length = fs_math::det::sqrt(0.5 * 0.5 + 0.2 * 0.2);
        for section in &mut laws[0] {
            section.resistance = MaterialFlowResistance::try_new(
                A * 0.5,
                0.0,
                Area::new(PIPE_AREA * 0.5 / length),
                Length::new(0.03),
                (0.0, 2300.0),
            )
            .unwrap();
        }
        let owner = LiquidGraphSources::try_new(
            cx,
            &law,
            &caloric,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [0.0, 0.0, -9.81],
            },
            caps(),
        )
        .unwrap();
        let problem = QuasiSteadyGraphProblem {
            topology: &graph,
            node_cells: inventory.node_cells(),
            reference_density: Density::new(1000.0),
            maximum_connections: 32,
        };
        let controls = QuasiSteadyGraphControls {
            solver: NewtonKrylovConfig {
                absolute_tolerance: 1e-9,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            maximum_cells: 16,
            maximum_newton_iterations: 32,
            maximum_endpoint_evaluations: 2048,
            energy_scale_j: 1.0,
            finite_difference_step: 1e-4,
            cell_tolerance_j: 1e-6,
            total_tolerance_j: 1e-6,
            maximum_neglected_power_w: 100.0,
            maximum_neglected_energy_j: 100.0,
        };
        let links = [GraphConvectiveLink {
            link: 9,
            fluid_cell: CellId(2),
            characteristic_length: Length::new(0.03),
            property_temperature: ConvectionPropertyTemperature::FilmMean,
            law: GraphConvectionLaw::Forced {
                branch: 0,
                correlation: CorrelationId::CircularDuctLaminarCwt,
                settings: developed_settings(),
            },
        }];
        let duration = 0.1;
        let (next, report) = advance_quasi_steady_graph_distributed(
            cx,
            &problem,
            &inventory,
            &initial,
            duration,
            controls,
            &[0.0; 7],
            |cx, candidate| {
                // These are source probes only. The native solved rates below
                // must independently pass the original resistance/heat cards.
                let resolved = owner
                    .resolve(cx, candidate, &[0.0; 3])
                    .map_err(|e| e.to_string())?;
                let hydraulics = resolved
                    .graph
                    .solve(
                        &[FixedPressure {
                            node: 0,
                            pressure: Pressure::new(0.0),
                        }],
                        hydraulic_controls(),
                        cx,
                    )
                    .map_err(|e| e.to_string())?;
                let heat = owner
                    .evaluate_quasi_steady_links(
                        cx,
                        candidate,
                        &resolved.sources,
                        &hydraulics,
                        duration,
                        &links,
                        16,
                    )
                    .map_err(|e| e.to_string())?;
                let mut conductances = BTreeMap::new();
                heat.convection
                    .apply_to(&mut conductances, || cx.checkpoint().is_err())
                    .map_err(|e| e.to_string())?;
                Ok::<_, String>(QuasiSteadyGraphResponse {
                    hydraulics,
                    mass_boundaries: vec![],
                    heat: StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: conductances,
                        payload: (resolved.sources, heat),
                    }
                    .into(),
                })
            },
        )
        .unwrap();
        let (sources, heat) = &report.heat;
        let sample = &heat.convection.samples()[0];
        let q = report.hydraulics.branches[0].flow.value();
        assert!(q > 0.0);
        assert!(heat.admission.actual_losses[0][0].reynolds > 0.0);
        assert!(heat.admission.actual_losses[0][0].reynolds < 2300.0);
        near(heat.admission.actual_losses[0][0].volume_flow_m3_s, q);
        assert!(next.temperature_k(CellId(2)).unwrap() > 300.0);
        assert!(next.temperature_k(CellId(6)).unwrap() < 400.0);
        near(next.time_s(), duration);
        assert_eq!(sample.link, 9);
        near(
            sample.property_temperature_k,
            0.5 * (next.temperature_k(CellId(2)).unwrap() + next.temperature_k(CellId(6)).unwrap()),
        );
        near(sample.conductance_w_k, 0.002 * 3.66 * 0.5 / 0.03);
        assert_eq!(
            sample.material.material().card_identity(),
            law.reference_state().material().card_identity()
        );
        let bulk = &sources.sections()[0][0].momentum.medium;
        near(
            bulk.density_kg_m3(),
            1000.0 - 0.5 * (next.temperature_k(CellId(2)).unwrap() - 300.0),
        );
        near(
            heat.admission.neglected_energy_j,
            duration * heat.admission.neglected_power_w,
        );
        near(report.neglected_power_w, heat.admission.neglected_power_w);

        let gravity = report.hydraulics.gravity.as_ref().unwrap();
        let sections = gravity.branch_sections().unwrap();
        assert_eq!(sections[0].len(), 2);
        assert_eq!(sections[0][0].end_m, [0.2, 0.0, 0.5]);
        assert_eq!(sections[0][1].start_m, sections[0][0].end_m);
        inventory.admit_gravity(cx, gravity).unwrap();
        assert!(report.maximum_cell_residual_j <= controls.cell_tolerance_j);
        assert!(report.thermal.residual_j.abs() <= controls.total_tolerance_j);
        near(report.thermal.external_heat_j, 0.0);
        let total_change: f64 = next
            .cells()
            .iter()
            .map(|(id, cell)| {
                cell.mass_kg * (cell.specific_h_j_kg - initial.cells()[id].specific_h_j_kg)
            })
            .sum();
        assert!(total_change.abs() < 1e-6);
        let wall = &next.cells()[&CellId(6)];
        let wall_change =
            wall.mass_kg * (wall.specific_h_j_kg - initial.cells()[&CellId(6)].specific_h_j_kg);
        assert!((wall_change + duration * sample.heat_into_fluid_w).abs() < 1e-6);
        near(initial.time_s(), 0.0);
        near(initial.temperature_k(CellId(2)).unwrap(), 300.0);
        near(initial.temperature_k(CellId(6)).unwrap(), 400.0);
    });
}

#[test]
fn actual_quasi_steady_reynolds_is_checked_after_an_admissible_source_probe() {
    context(|cx| {
        let card = card(false, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = thermal(&caloric, [300.0; 6], None, false);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(100.0);
        let owner = LiquidGraphSources::try_new(
            cx,
            &law,
            &caloric,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [0.0; 3],
            },
            caps(),
        )
        .unwrap();
        let resolved = owner.resolve(cx, &initial, &[0.0; 3]).unwrap();
        let flow = fixed_solution(cx, &resolved.graph, [1.0, 0.0]);
        near(flow.branches[0].flow.value(), 1e-5);
        match owner.admit_quasi_steady(cx, &initial, &resolved.sources, &flow, 0.1) {
            Err(SourceGraphError::Resistance {
                branch: 0,
                section: 0,
                source:
                    fs_airflow::graph::thermal::medium::series::hydraulic::ResistanceError::Reynolds {
                        value,
                        admitted,
                    },
            }) => {
                near(value, 150.0);
                assert_eq!(admitted, (0.0, 100.0));
            }
            other => panic!("actual quasi-steady Reynolds support was not enforced: {other:?}"),
        }
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[test]
fn quasi_steady_local_work_duration_and_native_interior_gravity_path_are_admitted() {
    context(|cx| {
        let card = buoyant_card();
        let (law, caloric) = source_pair(&card);
        let initial = thermal(
            &caloric,
            [300.0, 300.0, 280.0, 320.0, 300.0, 300.0],
            None,
            false,
        );
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(2300.0);
        let limits = SourceGraphLimits {
            maximum_neglected_section_energy_j: 0.001,
            ..caps()
        };
        let owner = LiquidGraphSources::try_new(
            cx,
            &law,
            &caloric,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [0.0, 0.0, -9.81],
            },
            limits,
        )
        .unwrap();
        let resolved = owner.resolve(cx, &initial, &[0.0; 3]).unwrap();
        // Ten pascals beyond the uniform-density hydrostatic frame drives all
        // three branches. The two local anomalies in branch zero cancel.
        let pressures = [10.0, -9810.0];
        let flow = fixed_solution(cx, &resolved.graph, pressures);
        let admitted = owner
            .admit_quasi_steady(cx, &initial, &resolved.sources, &flow, 0.01)
            .unwrap();
        near(flow.gravity.as_ref().unwrap().departure_heads_pa()[0], 0.0);
        near(admitted.neglected_buoyancy_power_w, 0.00981);
        near(admitted.neglected_hydraulic_power_w, 0.003);
        near(admitted.neglected_energy_j, 0.0001281);
        match owner.admit_quasi_steady(cx, &initial, &resolved.sources, &flow, 0.1) {
            Err(SourceGraphError::NeglectedWork { power_w, energy_j }) => {
                near(power_w, 0.01281);
                near(energy_j, 0.001281);
            }
            other => panic!("actual interval local work was not bounded: {other:?}"),
        }
        for duration in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                owner.admit_quasi_steady(cx, &initial, &resolved.sources, &flow, duration),
                Err(SourceGraphError::Invalid(_))
            ));
        }

        let gravity = resolved.graph.gravity().unwrap();
        let mut paths = gravity.branch_sections().unwrap().to_vec();
        paths[0][0].end_m[0] = 0.25;
        paths[0][1].start_m[0] = 0.25;
        let changed_graph = resolved
            .graph
            .clone()
            .with_distributed_boussinesq_gravity(
                law.reference_state(),
                &paths,
                gravity.positions_m(),
                gravity.acceleration_m_s2(),
                gravity.maximum_relative_density_departure(),
                gravity.maximum_sections().unwrap(),
                cx,
            )
            .unwrap();
        let changed = fixed_solution(cx, &changed_graph, pressures);
        // Node positions and every integrated head are identical, but the
        // solved branch followed a different interior path from stored fluid.
        assert_eq!(
            changed.gravity.as_ref().unwrap().pressure_heads_pa(),
            gravity.pressure_heads_pa()
        );
        assert!(matches!(
            owner.admit_quasi_steady(cx, &initial, &resolved.sources, &changed, 0.01),
            Err(SourceGraphError::Inventory(_))
        ));
        assert_eq!(initial.time_s(), 0.0);
    });
}
