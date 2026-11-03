//! Original developing-duct heat on the native material/flow/enthalpy owners.
//! Synthetic sources test implementation and do not establish calibration.

use super::*;
use fs_airflow::graph::thermal::transport::transient::material::convection::developing::{
    DevelopingWall, GraphDevelopingConvectionError, GraphDevelopingSection,
};
use fs_conduction::transient::deposition::{
    Reservoir, backward_euler::{BoundaryConductances, BoundaryStepCoefficients},
};
use fs_material::fluid_transport::FluidTransport;

fn diameter() -> f64 { (4.0 * PIPE_AREA / core::f64::consts::PI).sqrt() }
fn area() -> f64 { core::f64::consts::PI * diameter() * 0.5 }

fn section_fixture(caloric: &ResolvedHeatCapacityCurve, wall_k: f64) -> DepositionSimulation<'_> {
    let original = thermal(caloric, [300.0; 6], Some(wall_k), false);
    let mut schedule: Vec<_> = original.cells().values().cloned().map(|mut cell| {
        if cell.id == CellId(3) {
            cell.boundary.reservoir = Some(Reservoir {
                area_m2: area(), conductance_w_k: 0.0, temperature_k: wall_k,
            });
        }
        TimedEvent { time_s: 0.0, event: Event::Birth(cell) }
    }).collect();
    schedule.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
        area_m2: area(), ..original.links()[&9].clone()
    }) });
    DepositionSimulation::new(
        schedule, &[(MaterialId(1), caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 16, max_links: 16, max_events: 16, max_substeps: 16 },
    ).unwrap().advance_to(
        0.0, StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-8 }, || Ok(()),
    ).unwrap().0
}

fn sections() -> [GraphDevelopingSection; 2] {
    [
        GraphDevelopingSection {
            fluid_cell: CellId(2), wall: DevelopingWall::Link(9),
            diameter: Length::new(diameter()), forward_start: Length::new(0.0),
            reverse_start: Length::new(0.5),
            property_temperature: ConvectionPropertyTemperature::FilmMean,
        },
        GraphDevelopingSection {
            fluid_cell: CellId(3), wall: DevelopingWall::Reservoir,
            diameter: Length::new(diameter()), forward_start: Length::new(0.5),
            reverse_start: Length::new(0.0),
            property_temperature: ConvectionPropertyTemperature::FilmMean,
        },
    ]
}

fn hydraulic_controls() -> GraphSolveConfig {
    GraphSolveConfig {
        max_sweeps: 256, max_node_iterations: 96,
        absolute_flow_tolerance: VolumetricFlowRate::new(1e-14),
        relative_flow_tolerance: 1e-12,
    }
}

#[test]
fn g1_g3_developing_wall_sections_preserve_total_transfer_and_reverse_the_entrance() {
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let graph = graph();
        let laws = laws(2300.0);
        for wall_k in [280.0, 400.0] {
            let initial = section_fixture(&caloric, wall_k);
            let inventory = bind_inventory(cx, &graph, &initial);
            let owner = LiquidGraphSources::try_new(
                cx, &law, &caloric, &initial,
                SourceGraphBinding {
                    graph: &graph, inventory: &inventory, laws: &laws,
                    acceleration_m_s2: [0.0; 3],
                }, caps(),
            ).unwrap();
            for sign in [-1.0, 1.0] {
                let resolved = owner.resolve(cx, &initial, &[sign * 1e-5, sign * 1e-5, -sign * 1e-5]).unwrap();
                let hydraulics = resolved.graph.solve(&[
                    FixedPressure { node: 0, pressure: Pressure::new(sign) },
                    FixedPressure { node: 1, pressure: Pressure::new(0.0) },
                ], hydraulic_controls(), cx).unwrap();
                near(hydraulics.branches[0].flow.value(), sign * 1e-5);
                let bindings = sections();
                let heat = owner.evaluate_quasi_steady_developing_sections(
                    cx, &initial, &resolved.sources, &hydraulics, 0.25, &bindings, 2,
                ).unwrap();
                let samples = heat.convection.samples();
                let film = 0.5 * (300.0 + wall_k);
                let transport = law.at_temperature(Temperature::new(film)).unwrap().transport_state().unwrap();
                let re = transport.density_kg_m3() * 0.01 * diameter()
                    / transport.dynamic_viscosity_pa_s();
                let gz = re * transport.prandtl() * diameter(); // Whole heated length = 1 m.
                let nu = 3.66 + 0.0668 * gz / (1.0 + 0.04 * gz.powf(2.0 / 3.0));
                let total = core::f64::consts::PI * nu * transport.thermal_conductivity_w_m_k();
                near(samples.iter().map(|s| s.conductance_w_k).sum(), total);
                near(samples.iter().map(|s| s.heat_into_fluid_w).sum(), total * (wall_k - 300.0));
                assert_eq!(samples[0].conductance_w_k > samples[1].conductance_w_k, sign > 0.0);
                for sample in samples {
                    near(sample.property_temperature_k, film);
                    near(sample.area_m2, area());
                    near(sample.length_m, 0.5);
                    assert_eq!(sample.material.transport_state().unwrap(), transport);
                }
                let mut links = BTreeMap::new();
                let mut boundaries = BoundaryConductances::default();
                boundaries.ambient_w_k.insert(CellId(0), 7.0);
                heat.convection.apply_to(&mut links, &mut boundaries, || false).unwrap();
                near(links[&9], samples[0].conductance_w_k);
                near(boundaries.reservoir_w_k[&CellId(3)], samples[1].conductance_w_k);
                near(boundaries.ambient_w_k[&CellId(0)], 7.0);
                let prior = (links.clone(), boundaries.clone());
                assert!(heat.convection.apply_to(&mut links, &mut boundaries, || false).is_err());
                assert_eq!((links, boundaries), prior);
                let mut empty_links = BTreeMap::new();
                let mut conflict = BoundaryConductances::default();
                conflict.reservoir_w_k.insert(CellId(3), 1.0);
                assert!(heat.convection.apply_to(&mut empty_links, &mut conflict, || false).is_err());
                assert!(empty_links.is_empty());
                let mut empty = BoundaryConductances::default();
                assert!(heat.convection.apply_to(&mut empty_links, &mut empty, || true).is_err());
                assert!(empty_links.is_empty() && empty.reservoir_w_k.is_empty());
                let mut wrong_geometry = bindings;
                wrong_geometry[0].diameter = Length::new(1.01 * diameter());
                assert!(owner.evaluate_quasi_steady_developing_sections(
                    cx, &initial, &resolved.sources, &hydraulics, 0.25, &wrong_geometry, 2,
                ).is_err());
            }
            // Both nonzero endpoints and the average fit Hausen; crossing Re=0
            // still refuses, retaining the original forced-convection support.
            let average = [5e-6, 0.0, 5e-6];
            let resolved = owner.resolve(cx, &initial, &average).unwrap();
            let crossing = native_report(
                cx, &resolved.graph, &inventory, [-1e-5, 0.0, -1e-5],
                [300.5, 0.0, 300.5], average, 0.1,
            );
            assert!(matches!(owner.evaluate_developing_sections(
                cx, &initial, &resolved.sources, &crossing, &sections(), 2,
            ), Err(GraphDevelopingConvectionError::Convection(_))));
        }
    });
}

#[test]
fn g1_developing_links_and_reservoirs_share_native_momentum_advection_and_energy() {
    use fs_airflow::graph::thermal::transport::transient::inertial::{
        InertialGraphHeat, InertialThermalControls, InertialThermalProblem,
        ResolvedInertialGraph, advance_inertial_graph_distributed,
    };
    use fs_solver::NewtonKrylovConfig;
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = section_fixture(&caloric, 400.0);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(2300.0);
        let owner = LiquidGraphSources::try_new(
            cx, &law, &caloric, &initial,
            SourceGraphBinding {
                graph: &graph, inventory: &inventory, laws: &laws,
                acceleration_m_s2: [0.0; 3],
            }, caps(),
        ).unwrap();
        let initial_flows = [1e-5, 0.0, 1e-5];
        let first = owner.resolve(cx, &initial, &initial_flows).unwrap();
        let branches: Vec<_> = inventory.inertial_branches().iter().map(|b| b.unwrap()).collect();
        let boundaries = vec![FixedPressure { node: 0, pressure: Pressure::new(0.0) }];
        let momentum = InertialGraph::try_new(
            cx, &first.graph, Density::new(1000.0), &branches, &boundaries, &initial_flows, 0.0, 32,
        ).unwrap();
        let problem = InertialThermalProblem {
            momentum: &momentum, node_cells: inventory.node_cells(), maximum_connections: 32,
        };
        let controls = InertialThermalControls {
            solver: NewtonKrylovConfig {
                absolute_tolerance: 1e-9, relative_tolerance: 0.0, ..NewtonKrylovConfig::default()
            },
            maximum_cells: 16, maximum_newton_iterations: 32, maximum_endpoint_evaluations: 2048,
            energy_scale_j: 1.0, flow_scale_m3_s: 1e-5, finite_difference_step: 1e-4,
            cell_tolerance_j: 1e-6, total_thermal_tolerance_j: 1e-6,
            impulse_tolerance_pa_s: 1e-6, mechanical_tolerance_j: 1e-8,
            average_flow_tolerance_m3_s: 1e-12,
            momentum: InertialGraphControls {
                maximum_drag_step: 2.0, maximum_neglected_power_w: 100.0,
                maximum_neglected_energy_j: 100.0,
            },
        };
        let dt = 0.1;
        let (next_momentum, next, report) = advance_inertial_graph_distributed(
            cx, &problem, &inventory, &initial, dt, controls, &[0.0; 7], &initial_flows,
            |cx, candidate, flows| {
                let resolved = owner.resolve(cx, candidate, flows).map_err(|e| e.to_string())?;
                Ok::<_, String>(ResolvedInertialGraph {
                    graph: resolved.graph, pressure_boundaries: boundaries.clone(),
                    driving_pressure_pa: vec![1.0, 0.0, 1.0], payload: resolved.sources,
                })
            },
            |cx, candidate, sources, hydraulics| {
                let heat = owner.evaluate_developing_sections(
                    cx, candidate, sources, hydraulics, &sections(), 2,
                ).map_err(|e| e.to_string())?;
                let mut links = BTreeMap::new();
                let mut boundaries = BoundaryConductances::default();
                heat.convection.apply_to(&mut links, &mut boundaries, || cx.checkpoint().is_err())
                    .map_err(|e| e.to_string())?;
                Ok::<_, String>(InertialGraphHeat {
                    mass_boundaries: vec![],
                    heat: BoundaryStepCoefficients {
                        bulk: StepCoefficients {
                            sources_w: BTreeMap::new(), conductances_w_k: links, payload: heat,
                        }, boundaries,
                    },
                })
            },
        ).unwrap();
        near(next.time_s(), dt);
        near(next_momentum.accepted_time_s(), dt);
        assert!(next.temperature_k(CellId(2)).unwrap() > 300.0);
        assert!(next.temperature_k(CellId(3)).unwrap() > 300.0);
        assert!(next.temperature_k(CellId(6)).unwrap() < 400.0);
        let samples = report.heat.convection.samples();
        // Only the prescribed reservoir is external; native finite-wall and
        // enthalpy advection transfers cancel inside this closed fluid loop.
        let external = dt * samples[1].heat_into_fluid_w;
        let stored: f64 = next.cells().iter().map(|(id, cell)| {
            cell.mass_kg * (cell.specific_h_j_kg - initial.cells()[id].specific_h_j_kg)
        }).sum();
        assert!((stored - external).abs() < 2e-6);
        assert!((report.thermal.external_heat_j - external).abs() < 2e-6);
        for sample in samples {
            near(sample.fluid_temperature_k, next.temperature_k(sample.fluid_cell).unwrap());
            near(sample.property_temperature_k,
                0.5 * (sample.wall_temperature_k + sample.fluid_temperature_k));
        }
        assert!(report.maximum_cell_residual_j <= controls.cell_tolerance_j);
        assert!(report.hydraulics.mechanical_balance_residual_j.abs() <= controls.mechanical_tolerance_j);
        near(initial.time_s(), 0.0);
        near(momentum.accepted_time_s(), 0.0);
    });
}
