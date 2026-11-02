//! Actual liquid section loss enters native finite-cell heat exactly once.
//! Synthetic constitutive data test coupling and conservation, not calibration.

use super::*;
use fs_airflow::graph::thermal::transport::transient::{
    inertial::{
        InertialGraphHeat, InertialThermalControls, InertialThermalProblem,
        ResolvedInertialGraph, advance_inertial_graph_distributed,
    },
    material::heating::GraphHeatingError,
};
use fs_solver::NewtonKrylovConfig;

#[test]
fn g1_g3_actual_section_viscosity_sets_signed_native_loss_and_local_heat() {
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = thermal(
            &caloric,
            [300.0, 300.0, 300.0, 400.0, 300.0, 300.0],
            None,
            false,
        );
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(1e6);
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
        let heating = owner.dissipative_heating(cx).unwrap();
        for sign in [-1.0, 0.0, 1.0] {
            let resolved = owner
                .resolve(cx, &initial, &[sign * 0.001, sign * 0.00075, -sign * 0.00075])
                .unwrap();
            let hydraulics = resolved
                .graph
                .solve(
                    &[
                        FixedPressure {
                            node: 0,
                            pressure: Pressure::new(sign * 75.0),
                        },
                        FixedPressure {
                            node: 1,
                            pressure: Pressure::new(0.0),
                        },
                    ],
                    GraphSolveConfig {
                        max_sweeps: 128,
                        max_node_iterations: 80,
                        absolute_flow_tolerance: VolumetricFlowRate::new(1e-12),
                        relative_flow_tolerance: 1e-10,
                    },
                    cx,
                )
                .unwrap();
            near(hydraulics.branches[0].flow.value(), sign * 0.001);
            near(hydraulics.branches[1].flow.value(), sign * 0.00075);
            near(hydraulics.branches[2].flow.value(), -sign * 0.00075);
            let heat = heating
                .evaluate_quasi_steady(cx, &initial, &resolved.sources, &hydraulics, 0.5)
                .unwrap();
            // Independently: R = 50 kPa s/m3 and 25 kPa s/m3 in the
            // unequal-temperature halves, with 100 kPa s/m3 in each other pipe.
            let powers = if sign == 0.0 {
                [0.0; 4]
            } else {
                [0.05, 0.025, 0.05625, 0.05625]
            };
            assert_eq!(heat.samples().len(), 4);
            for (sample, ((cell, branch, section), power)) in heat.samples().iter().zip(
                [(2, 0, 0), (3, 0, 1), (4, 1, 0), (5, 2, 0)]
                    .into_iter()
                    .zip(powers),
            ) {
                assert_eq!((sample.cell, sample.branch, sample.section), (CellId(cell), branch, section));
                near(sample.power_w, power);
                near(sample.heat_j, 0.5 * power);
                assert!(sample.power_w >= 0.0);
                if sign != 0.0 {
                    near(heat.sources_w()[&sample.cell], power);
                }
            }
            near(heat.duration_s(), 0.5);
            near(heat.thermal_heat_j(), 0.5 * powers.iter().sum::<f64>());
            near(heat.mechanical_loss_j(), heat.thermal_heat_j());
            assert!(heat.transfer_residual_j().abs() < 1e-12);
            assert!(std::ptr::eq(
                heat.source_sections(),
                resolved.sources.sections(),
            ));
            assert!(!heat.sources_w().contains_key(&CellId(0)));
            assert!(!heat.sources_w().contains_key(&CellId(1)));
            if sign == 0.0 {
                assert!(heat.sources_w().is_empty());
            } else {
                // Additional heat can coexist with another cell's source.
                // Applying the same owner twice must leave the map unchanged.
                let mut additional = BTreeMap::from([(CellId(0), 7.0)]);
                heat.apply_to(cx, &mut additional).unwrap();
                near(additional[&CellId(0)], 7.0);
                let before = additional.clone();
                assert!(matches!(
                    heat.apply_to(cx, &mut additional),
                    Err(GraphHeatingError::CellSourceAlreadyOwned(_))
                ));
                assert_eq!(additional, before);
            }
            let changed = thermal(
                &caloric,
                [300.0, 300.0, 301.0, 400.0, 300.0, 300.0],
                None,
                false,
            );
            assert!(matches!(
                heating.evaluate_quasi_steady(
                    cx,
                    &changed,
                    &resolved.sources,
                    &hydraulics,
                    0.5,
                ),
                Err(GraphHeatingError::Source(SourceGraphError::Invalid(_)))
            ));
        }
    });
}

#[test]
fn g1_native_unforced_loop_converts_kinetic_loss_to_heat_once_at_actual_average_flow() {
    context(|cx| {
        for sign in [-1.0, 1.0] {
            let card = card(true, 1.0);
            let (law, caloric) = source_pair(&card);
            let initial = thermal(&caloric, [300.0; 6], None, false);
            let graph = graph();
            let inventory = bind_inventory(cx, &graph, &initial);
            let laws = laws(1e6);
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
            let heating = owner.dissipative_heating(cx).unwrap();
            let initial_flows = [sign * 0.01, 0.0, sign * 0.01];
            let first = owner.resolve(cx, &initial, &initial_flows).unwrap();
            let branches: Vec<_> = inventory
                .inertial_branches()
                .iter()
                .map(|branch| branch.unwrap())
                .collect();
            let boundaries = vec![FixedPressure {
                node: 0,
                pressure: Pressure::new(0.0),
            }];
            let momentum = InertialGraph::try_new(
                cx,
                &first.graph,
                Density::new(1000.0),
                &branches,
                &boundaries,
                &initial_flows,
                0.0,
                32,
            )
            .unwrap();
            // Two unit-length pipes each contain 1 kg moving at 10 m/s.
            near(momentum.kinetic_energy_j(), 100.0);
            let problem = InertialThermalProblem {
                momentum: &momentum,
                node_cells: inventory.node_cells(),
                maximum_connections: 32,
            };
            let controls = InertialThermalControls {
                solver: NewtonKrylovConfig {
                    absolute_tolerance: 1e-9,
                    relative_tolerance: 0.0,
                    ..NewtonKrylovConfig::default()
                },
                maximum_cells: 16,
                maximum_newton_iterations: 48,
                maximum_endpoint_evaluations: 4096,
                energy_scale_j: 1.0,
                flow_scale_m3_s: 0.01,
                finite_difference_step: 1e-4,
                cell_tolerance_j: 1e-7,
                total_thermal_tolerance_j: 1e-7,
                impulse_tolerance_pa_s: 1e-6,
                mechanical_tolerance_j: 1e-7,
                average_flow_tolerance_m3_s: 1e-11,
                momentum: InertialGraphControls {
                    maximum_drag_step: 2.0,
                    maximum_neglected_power_w: 100.0,
                    maximum_neglected_energy_j: 100.0,
                },
            };
            let (next_momentum, next, report) = advance_inertial_graph_distributed(
                cx,
                &problem,
                &inventory,
                &initial,
                0.25,
                controls,
                &[0.0; 6],
                &initial_flows,
                |cx, candidate, flows| {
                    let resolved = owner
                        .resolve(cx, candidate, flows)
                        .map_err(GraphHeatingError::Source)?;
                    Ok::<_, GraphHeatingError>(ResolvedInertialGraph {
                        graph: resolved.graph,
                        pressure_boundaries: boundaries.clone(),
                        driving_pressure_pa: vec![0.0; 3],
                        payload: resolved.sources,
                    })
                },
                |cx, candidate, sources, hydraulics| {
                    let heat = heating.evaluate_inertial(cx, candidate, sources, hydraulics)?;
                    let mut sources_w = BTreeMap::new();
                    heat.apply_to(cx, &mut sources_w)?;
                    Ok::<_, GraphHeatingError>(InertialGraphHeat {
                        mass_boundaries: vec![],
                        heat: StepCoefficients {
                            sources_w,
                            conductances_w_k: BTreeMap::new(),
                            payload: heat,
                        }
                        .into(),
                    })
                },
            )
            .unwrap();
            let stored: f64 = next
                .cells()
                .iter()
                .map(|(id, cell)| {
                    cell.mass_kg * (cell.specific_h_j_kg - initial.cells()[id].specific_h_j_kg)
                })
                .sum();
            let kinetic_change = next_momentum.kinetic_energy_j() - momentum.kinetic_energy_j();
            assert!(report.hydraulics.hydraulic_loss_j > 1.0);
            assert!(kinetic_change < 0.0);
            assert!((stored - report.hydraulics.hydraulic_loss_j).abs() <= 1e-7);
            assert!((stored + kinetic_change).abs() <= 2e-7);
            near(report.heat.mechanical_loss_j(), report.hydraulics.hydraulic_loss_j);
            near(report.heat.thermal_heat_j(), report.hydraulics.hydraulic_loss_j);
            near(report.thermal.external_heat_j, report.heat.thermal_heat_j());
            near(report.hydraulics.driving_work_j, 0.0);
            near(report.hydraulics.reduced_boundary_work_j, 0.0);
            near(report.hydraulics.anomalous_gravity_work_j, 0.0);
            assert!(
                report.heat.admission().neglected_hydraulic_energy_j
                    > report.heat.thermal_heat_j(),
                "the endpoint speed envelope is an admission bound, never heat",
            );
            assert!(next.temperature_k(CellId(2)).unwrap() > 300.0);
            let actual_t = next.temperature_k(CellId(2)).unwrap();
            let source = &report.sources.sections()[0][0].momentum;
            near(source.medium.dynamic_viscosity_pa_s(), 0.005 - 1e-5 * actual_t);
            assert!(source.medium.dynamic_viscosity_pa_s() < 0.002);
            near(
                source.evaluation.volume_flow_m3_s,
                report.hydraulics.average_volume_flows_m3_s[0],
            );
            assert!(std::ptr::eq(
                report.heat.source_sections(),
                report.sources.sections(),
            ));
            near(next.time_s(), 0.25);
            near(next_momentum.accepted_time_s(), 0.25);
            near(initial.time_s(), 0.0);
            near(momentum.accepted_time_s(), 0.0);
            near(initial.temperature_k(CellId(2)).unwrap(), 300.0);
        }
    });
}

#[test]
fn g4_original_density_span_and_native_mechanical_debit_cannot_be_substituted() {
    context(|cx| {
        // rho(300 K) is exactly 1000, yet its original source span is not
        // constant. Qualifying only the reference or current state would pass.
        let variable_card = buoyant_card();
        let (variable_law, variable_caloric) = source_pair(&variable_card);
        let variable_initial = thermal(&variable_caloric, [300.0; 6], None, false);
        let variable_graph = graph();
        let variable_inventory = bind_inventory(cx, &variable_graph, &variable_initial);
        let laws = laws(1e6);
        let variable_owner = LiquidGraphSources::try_new(
            cx,
            &variable_law,
            &variable_caloric,
            &variable_initial,
            SourceGraphBinding {
                graph: &variable_graph,
                inventory: &variable_inventory,
                laws: &laws,
                acceleration_m_s2: [0.0; 3],
            },
            caps(),
        )
        .unwrap();
        assert!(matches!(
            variable_owner.dissipative_heating(cx),
            Err(GraphHeatingError::Source(SourceGraphError::LiquidCaloric(_)))
        ));

        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = thermal(&caloric, [300.0; 6], None, false);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let binding = || SourceGraphBinding {
            graph: &graph,
            inventory: &inventory,
            laws: &laws,
            acceleration_m_s2: [0.0; 3],
        };
        let owner = LiquidGraphSources::try_new(
            cx, &law, &caloric, &initial, binding(), caps(),
        )
        .unwrap();
        let heating = owner.dissipative_heating(cx).unwrap();
        let flows = [0.001, 0.0, 0.001];
        let resolved = owner.resolve(cx, &initial, &flows).unwrap();
        let report = native_report(
            cx,
            &resolved.graph,
            &inventory,
            flows,
            [100.0, 0.0, 100.0],
            flows,
            0.25,
        );
        heating
            .evaluate_inertial(cx, &initial, &resolved.sources, &report)
            .unwrap();
        let mut substituted = report.clone();
        substituted.hydraulic_loss_j += 1.0;
        assert!(matches!(
            heating.evaluate_inertial(cx, &initial, &resolved.sources, &substituted),
            Err(GraphHeatingError::TransferBalance { .. })
        ));

        // This opt-in does not relax the original source applicability caps.
        let mut zero_caps = caps();
        zero_caps.maximum_neglected_section_power_w = 0.0;
        zero_caps.maximum_neglected_section_energy_j = 0.0;
        let bounded = LiquidGraphSources::try_new(
            cx, &law, &caloric, &initial, binding(), zero_caps,
        )
        .unwrap();
        let bounded_sources = bounded.resolve(cx, &initial, &flows).unwrap();
        assert!(matches!(
            bounded.dissipative_heating(cx).unwrap().evaluate_inertial(
                cx,
                &initial,
                &bounded_sources.sources,
                &report,
            ),
            Err(GraphHeatingError::Source(SourceGraphError::NeglectedWork { .. }))
        ));
    });
}
