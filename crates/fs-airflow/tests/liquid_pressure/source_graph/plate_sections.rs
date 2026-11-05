//! Physical plate strips on the original flow, material and finite heat owners.
//! Synthetic source data verify numerical behavior, not experimental accuracy.

use super::*;
use fs_airflow::graph::thermal::transport::transient::material::convection::{
    GraphConvectionError, PlateFlowDirection,
    reservoir::{GraphConvectiveReservoir, GraphReservoirConvectionError},
};
use fs_conduction::transient::deposition::{
    Reservoir,
    backward_euler::{BoundaryConductances, BoundaryStepCoefficients},
};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

const STRIP_AREA: f64 = 0.02;

fn fixture(caloric: &ResolvedHeatCapacityCurve) -> DepositionSimulation<'_> {
    let original = thermal(caloric, [300.0; 6], Some(400.0), false);
    let mut schedule: Vec<_> = original
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            if cell.id == CellId(3) {
                cell.boundary.reservoir = Some(Reservoir {
                    area_m2: STRIP_AREA,
                    conductance_w_k: 0.0,
                    temperature_k: 400.0,
                });
            }
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    schedule.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            area_m2: STRIP_AREA,
            ..original.links()[&9].clone()
        }),
    });
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 16,
            max_links: 16,
            max_events: 16,
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

fn bindings(direction: PlateFlowDirection) -> (GraphConvectiveLink, GraphConvectiveReservoir) {
    let starts = match direction {
        PlateFlowDirection::AlongBranch => [0.0, 0.5],
        PlateFlowDirection::AgainstBranch => [0.5, 0.0],
    };
    let law = |start| GraphConvectionLaw::ForcedPlateSection {
        branch: 0,
        correlation: CorrelationId::FlatPlateLaminarAverage,
        start: Length::new(start),
        direction,
    };
    (
        GraphConvectiveLink {
            link: 9,
            fluid_cell: CellId(2),
            characteristic_length: Length::new(starts[0] + 0.5),
            property_temperature: ConvectionPropertyTemperature::FilmMean,
            law: law(starts[0]),
        },
        GraphConvectiveReservoir {
            fluid_cell: CellId(3),
            characteristic_length: Length::new(starts[1] + 0.5),
            property_temperature: ConvectionPropertyTemperature::FilmMean,
            law: law(starts[1]),
        },
    )
}

// Independent laminar cumulative-transfer algebra. For h_bar(x)=C/sqrt(x),
// a finite strip has h=C/(sqrt(start)+sqrt(end)); no production evaluator is used.
fn laminar_h(transport: FluidTransportState, speed: f64, start: f64, end: f64) -> f64 {
    let c = 0.664
        * transport.thermal_conductivity_w_m_k()
        * (transport.density_kg_m3() * speed / transport.dynamic_viscosity_pa_s()).sqrt()
        * transport.prandtl().powf(1.0 / 3.0);
    c / (start.sqrt() + end.sqrt())
}

#[test]
fn g1_g3_native_plate_strips_preserve_total_heat_and_rebind_the_physical_leading_edge() {
    context(|cx| {
        let card = card(false, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = fixture(&caloric);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(2300.0);
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
        let q = 1e-5;
        for (sign, direction) in [
            (1.0, PlateFlowDirection::AlongBranch),
            (-1.0, PlateFlowDirection::AgainstBranch),
        ] {
            let resolved = owner
                .resolve(cx, &initial, &[sign * q, sign * q, -sign * q])
                .unwrap();
            let hydraulics = resolved
                .graph
                .solve(
                    &[
                        FixedPressure {
                            node: 0,
                            pressure: Pressure::new(sign),
                        },
                        FixedPressure {
                            node: 1,
                            pressure: Pressure::new(0.0),
                        },
                    ],
                    GraphSolveConfig {
                        max_sweeps: 256,
                        max_node_iterations: 96,
                        absolute_flow_tolerance: VolumetricFlowRate::new(1e-14),
                        relative_flow_tolerance: 1e-12,
                    },
                    cx,
                )
                .unwrap();
            near(hydraulics.branches[0].flow.value(), sign * q);
            let (link, reservoir) = bindings(direction);
            let links = owner
                .evaluate_quasi_steady_links(
                    cx,
                    &initial,
                    &resolved.sources,
                    &hydraulics,
                    0.1,
                    &[link],
                    1,
                )
                .unwrap();
            let reservoirs = owner
                .evaluate_quasi_steady_reservoirs(
                    cx,
                    &initial,
                    &resolved.sources,
                    &hydraulics,
                    0.1,
                    &[reservoir],
                    1,
                )
                .unwrap();
            let first = &links.convection.samples()[0];
            let second = &reservoirs.convection.samples()[0];
            for (sample, transport, end) in [
                (
                    first.evaluation.as_ref().unwrap(),
                    first.thermal.transport(),
                    link.characteristic_length,
                ),
                (
                    second.evaluation.as_ref().unwrap(),
                    second.thermal.transport(),
                    reservoir.characteristic_length,
                ),
            ] {
                let start = end.value() - 0.5;
                near(
                    sample.coefficient.value.value(),
                    laminar_h(transport, q / PIPE_AREA, start, end.value()),
                );
                let section = sample.section.as_ref().unwrap();
                near(section.start.value(), start);
                near(section.end.value(), end.value());
                assert_eq!(section.upstream.is_some(), start > 0.0);
                if start > 0.0 {
                    assert!(sample.effective_nusselt() < sample.nusselt.evidence().value);
                }
            }
            let whole_h = laminar_h(first.thermal.transport(), q / PIPE_AREA, 0.0, 1.0);
            near(
                first.conductance_w_k + second.conductance_w_k,
                2.0 * STRIP_AREA * whole_h,
            );
            near(
                first.heat_into_fluid_w + second.heat_into_fluid_w,
                100.0 * 2.0 * STRIP_AREA * whole_h,
            );
            near(first.property_temperature_k, 350.0);
            near(second.property_temperature_k, 350.0);
            assert_eq!(
                first.material.material().card_identity(),
                law.reference_state().material().card_identity()
            );
            assert_eq!(
                second.material.material().card_identity(),
                first.material.material().card_identity()
            );
        }
    });
}

#[test]
fn g4_native_plate_sections_refuse_reversal_and_hidden_upstream_regime_failure() {
    context(|cx| {
        let card = card(false, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = fixture(&caloric);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(2300.0);
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
        let q = 1e-5;
        let resolved = owner.resolve(cx, &initial, &[q, 0.0, q]).unwrap();
        let (link, reservoir) = bindings(PlateFlowDirection::AlongBranch);
        // All source/drag ranges admit these trajectories. The mean is positive,
        // but the first trajectory crosses the physical leading edge; the second
        // has Re_end=1.5 and Re_start=0.75 at its old endpoint.
        for old in [-q, 3e-9] {
            let new = 2.0 * q - old;
            let drive = 1e6 * (new - old) / 0.1 + 1e5 * q;
            let report = native_report(
                cx,
                &resolved.graph,
                &inventory,
                [old, 0.0, old],
                [drive, 0.0, drive],
                [q, 0.0, q],
                0.1,
            );
            owner
                .admit_actual(cx, &initial, &resolved.sources, &report)
                .unwrap();
            let result = owner.evaluate_reservoirs(
                cx,
                &initial,
                &resolved.sources,
                &report,
                &[reservoir],
                1,
            );
            if old < 0.0 {
                assert!(matches!(
                    result,
                    Err(GraphReservoirConvectionError::Source(
                        SourceGraphError::Invalid(_)
                    ))
                ));
                assert!(matches!(
                    owner.evaluate_links(cx, &initial, &resolved.sources, &report, &[link], 1,),
                    Err(GraphConvectionError::Source(SourceGraphError::Invalid(_)))
                ));
            } else {
                assert!(matches!(
                    result,
                    Err(GraphReservoirConvectionError::Envelope {
                        cell: CellId(3),
                        ..
                    })
                ));
            }
        }
        let steady = native_report(
            cx,
            &resolved.graph,
            &inventory,
            [q, 0.0, q],
            [1.0, 0.0, 1.0],
            [q, 0.0, q],
            0.1,
        );
        let mut wrong_film = link;
        wrong_film.property_temperature = ConvectionPropertyTemperature::Reservoir;
        assert!(matches!(
            owner.evaluate_links(cx, &initial, &resolved.sources, &steady, &[wrong_film], 1,),
            Err(GraphConvectionError::Convection(_))
        ));
        let mut wrong_section = reservoir;
        wrong_section.characteristic_length = Length::new(0.5); // start == end
        assert!(matches!(
            owner.evaluate_reservoirs(
                cx,
                &initial,
                &resolved.sources,
                &steady,
                &[wrong_section],
                1,
            ),
            Err(GraphReservoirConvectionError::Convection(_))
        ));
        near(initial.time_s(), 0.0);
        near(initial.temperature_k(CellId(2)).unwrap(), 300.0);
        near(initial.temperature_k(CellId(6)).unwrap(), 400.0);
        near(initial.links()[&9].conductance_w_k, 0.0);
    });
}

#[test]
fn g1_native_plate_heat_updates_momentum_film_properties_and_endpoint_energy_together() {
    use fs_airflow::graph::thermal::transport::transient::inertial::{
        InertialGraphHeat, InertialThermalControls, InertialThermalProblem, ResolvedInertialGraph,
        advance_inertial_graph_distributed,
    };
    use fs_solver::NewtonKrylovConfig;
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = fixture(&caloric);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(2300.0);
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
        let initial_flows = [1e-5, 0.0, 1e-5];
        let first = owner.resolve(cx, &initial, &initial_flows).unwrap();
        let branches: Vec<_> = inventory
            .inertial_branches()
            .iter()
            .map(|b| b.unwrap())
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
            maximum_newton_iterations: 32,
            maximum_endpoint_evaluations: 2048,
            energy_scale_j: 1.0,
            flow_scale_m3_s: 1e-5,
            finite_difference_step: 1e-4,
            cell_tolerance_j: 1e-6,
            total_thermal_tolerance_j: 1e-6,
            impulse_tolerance_pa_s: 1e-6,
            mechanical_tolerance_j: 1e-8,
            average_flow_tolerance_m3_s: 1e-12,
            momentum: InertialGraphControls {
                maximum_drag_step: 2.0,
                maximum_neglected_power_w: 100.0,
                maximum_neglected_energy_j: 100.0,
            },
        };
        let (link, reservoir) = bindings(PlateFlowDirection::AlongBranch);
        let dt = 0.1;
        let (next_momentum, next, report) = advance_inertial_graph_distributed(
            cx,
            &problem,
            &inventory,
            &initial,
            dt,
            controls,
            &[0.0; 7],
            &initial_flows,
            |cx, candidate, flows| {
                let resolved = owner
                    .resolve(cx, candidate, flows)
                    .map_err(|e| e.to_string())?;
                Ok::<_, String>(ResolvedInertialGraph {
                    graph: resolved.graph,
                    pressure_boundaries: boundaries.clone(),
                    driving_pressure_pa: vec![1.0, 0.0, 1.0],
                    payload: resolved.sources,
                })
            },
            |cx, candidate, sources, hydraulics| {
                let links = owner
                    .evaluate_links(cx, candidate, sources, hydraulics, &[link], 1)
                    .map_err(|e| e.to_string())?;
                let reservoirs = owner
                    .evaluate_reservoirs(cx, candidate, sources, hydraulics, &[reservoir], 1)
                    .map_err(|e| e.to_string())?;
                let mut conductances = BTreeMap::new();
                let mut boundaries = BoundaryConductances::default();
                links
                    .convection
                    .apply_to(&mut conductances, || cx.checkpoint().is_err())
                    .map_err(|e| e.to_string())?;
                reservoirs
                    .convection
                    .apply_to(&mut boundaries, || cx.checkpoint().is_err())
                    .map_err(|e| e.to_string())?;
                Ok::<_, String>(InertialGraphHeat {
                    mass_boundaries: vec![],
                    heat: BoundaryStepCoefficients {
                        bulk: StepCoefficients {
                            sources_w: BTreeMap::new(),
                            conductances_w_k: conductances,
                            payload: (links, reservoirs),
                        },
                        boundaries,
                    },
                })
            },
        )
        .unwrap();
        let (links, reservoirs) = &report.heat;
        let wall = &links.convection.samples()[0];
        let bath = &reservoirs.convection.samples()[0];
        let q = report.hydraulics.average_volume_flows_m3_s[0];
        assert!(q > initial_flows[0]); // Heating reduces the original bulk drag.
        for (id, material, film, coefficient, start, end) in [
            (
                wall.fluid_cell,
                &wall.material,
                wall.property_temperature_k,
                wall.evaluation.as_ref().unwrap().coefficient.value.value(),
                0.0,
                0.5,
            ),
            (
                bath.fluid_cell,
                &bath.material,
                bath.property_temperature_k,
                bath.evaluation.as_ref().unwrap().coefficient.value.value(),
                0.5,
                1.0,
            ),
        ] {
            let transport = material.transport_state().unwrap();
            assert!(next.temperature_k(id).unwrap() > 300.0);
            near(transport.dynamic_viscosity_pa_s(), 0.005 - 1e-5 * film);
            near(coefficient, laminar_h(transport, q / PIPE_AREA, start, end));
        }
        near(
            wall.wall_temperature_k,
            next.temperature_k(CellId(6)).unwrap(),
        );
        near(
            wall.property_temperature_k,
            0.5 * (wall.wall_temperature_k + next.temperature_k(CellId(2)).unwrap()),
        );
        near(
            bath.property_temperature_k,
            0.5 * (400.0 + next.temperature_k(CellId(3)).unwrap()),
        );
        let frozen = law
            .at_temperature(Temperature::new(350.0))
            .unwrap()
            .transport_state()
            .unwrap();
        let frozen_h = laminar_h(frozen, initial_flows[0] / PIPE_AREA, 0.0, 0.5);
        assert!(
            (wall.evaluation.as_ref().unwrap().coefficient.value.value() - frozen_h).abs()
                > 1e-9 * frozen_h
        );
        let stored: f64 = next
            .cells()
            .iter()
            .map(|(id, cell)| {
                cell.mass_kg * (cell.specific_h_j_kg - initial.cells()[id].specific_h_j_kg)
            })
            .sum();
        let wall_change = next.cells()[&CellId(6)].mass_kg
            * (next.cells()[&CellId(6)].specific_h_j_kg
                - initial.cells()[&CellId(6)].specific_h_j_kg);
        let external = dt * bath.heat_into_fluid_w;
        assert!((stored - external).abs() < 2e-6);
        assert!((report.thermal.external_heat_j - external).abs() < 2e-6);
        assert!((wall_change + dt * wall.heat_into_fluid_w).abs() < 1e-6);
        assert!(report.maximum_cell_residual_j <= controls.cell_tolerance_j);
        assert!(
            report.hydraulics.mechanical_balance_residual_j.abs()
                <= controls.mechanical_tolerance_j
        );
        near(next.time_s(), dt);
        near(next_momentum.accepted_time_s(), dt);
        near(initial.time_s(), 0.0);
        near(momentum.accepted_time_s(), 0.0);
    });
}
