//! Prescribed wall reservoirs on the original liquid graph and thermal clock.
//! Synthetic source/correlation data test implementation, not calibration.

use super::*;
use fs_airflow::graph::thermal::medium::convection::NaturalConvectionGeometry;
use fs_airflow::graph::thermal::transport::transient::{
    coupled::{
        QuasiSteadyGraphControls, QuasiSteadyGraphProblem, QuasiSteadyGraphResponse,
        advance_quasi_steady_graph_distributed,
    },
    material::convection::reservoir::{GraphConvectiveReservoir, GraphReservoirConvectionError},
};
use fs_conduction::transient::deposition::{
    Reservoir,
    backward_euler::{BoundaryConductances, BoundaryStepCoefficients},
};
use fs_material::fluid_transport::FluidTransport;
use fs_solver::NewtonKrylovConfig;

fn with_reservoir(
    caloric: &ResolvedHeatCapacityCurve,
    temperatures: [f64; 6],
    reservoir: Option<Reservoir>,
) -> DepositionSimulation<'_> {
    let original = thermal(caloric, temperatures, None, false);
    let schedule = original.cells().values().cloned().map(|mut cell| {
        if cell.id == CellId(2) {
            cell.boundary.reservoir = reservoir;
        }
        TimedEvent { time_s: 0.0, event: Event::Birth(cell) }
    }).collect();
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 16, max_links: 16, max_events: 16, max_substeps: 16 },
    ).unwrap().advance_to(
        0.0,
        StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-8 },
        || Ok(()),
    ).unwrap().0
}

fn wall(temperature_k: f64) -> Reservoir {
    Reservoir { area_m2: 0.002, conductance_w_k: 0.0, temperature_k }
}

fn forced() -> GraphConvectiveReservoir {
    GraphConvectiveReservoir {
        fluid_cell: CellId(2),
        characteristic_length: Length::new(0.03),
        property_temperature: ConvectionPropertyTemperature::FilmMean,
        law: GraphConvectionLaw::Forced {
            branch: 0,
            correlation: CorrelationId::CircularDuctLaminarCwt,
            settings: developed_settings(),
        },
    }
}

fn hydraulic_controls() -> GraphSolveConfig {
    GraphSolveConfig {
        max_sweeps: 256,
        max_node_iterations: 96,
        absolute_flow_tolerance: VolumetricFlowRate::new(1e-14),
        relative_flow_tolerance: 1e-12,
    }
}

#[test]
fn g1_g3_prescribed_wall_keeps_full_heat_drop_and_original_film_source_under_reversal() {
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let graph = graph();
        let laws = laws(2300.0);
        for wall_temperature in [280.0, 400.0] {
            let initial = with_reservoir(&caloric, [300.0; 6], Some(wall(wall_temperature)));
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
                let flow = resolved.graph.solve(
                    &[
                        FixedPressure { node: 0, pressure: Pressure::new(sign) },
                        FixedPressure { node: 1, pressure: Pressure::new(0.0) },
                    ], hydraulic_controls(), cx,
                ).unwrap();
                near(flow.branches[0].flow.value(), sign * 1e-5);
                let heat = owner.evaluate_quasi_steady_reservoirs(
                    cx, &initial, &resolved.sources, &flow, 0.25, &[forced()], 1,
                ).unwrap();
                let sample = &heat.convection.samples()[0];
                let film = 0.5 * (300.0 + wall_temperature);
                // Fully developed circular CWT Nu=3.66, original k=0.5 W/mK.
                let conductance = 0.002 * 3.66 * 0.5 / 0.03;
                near(sample.wall_temperature_k, wall_temperature);
                near(sample.fluid_temperature_k, 300.0);
                near(sample.property_temperature_k, film);
                near(sample.thermal.temperature_k(), film);
                near(sample.area_m2, 0.002);
                near(sample.conductance_w_k, conductance);
                near(sample.heat_into_fluid_w, conductance * (wall_temperature - 300.0));
                near(sample.thermal.transport().dynamic_viscosity_pa_s(), 0.005 - 1e-5 * film);
                assert_eq!(
                    sample.material.transport_state().unwrap(),
                    law.at_temperature(Temperature::new(film)).unwrap().transport_state().unwrap(),
                );
                near(heat.convection.reservoir_conductances_w_k()[&CellId(2)], conductance);
                assert_eq!(heat.convection.reservoir_conductances_w_k().len(), 1);
                let mut native = BoundaryConductances::default();
                native.ambient_w_k.insert(CellId(0), 7.0);
                heat.convection.apply_to(&mut native, || false).unwrap();
                near(native.ambient_w_k[&CellId(0)], 7.0);
                near(native.reservoir_w_k[&CellId(2)], conductance);
                let accepted = native.clone();
                assert!(heat.convection.apply_to(&mut native, || false).is_err());
                assert_eq!(native, accepted);
                let mut cancelled = BoundaryConductances::default();
                assert!(heat.convection.apply_to(&mut cancelled, || true).is_err());
                assert_eq!(cancelled, BoundaryConductances::default());
                near(initial.cells()[&CellId(2)].boundary.reservoir.unwrap().conductance_w_k, 0.0);
            }
        }
    });
}

// Independent Churchill-Chu evaluation for the original liquid card:
// rho(T)=1000-0.5(T-300), beta=0.5/rho, mu=0.002, Cp=2000, k=0.5.
fn natural_conductance(fluid_temperature: f64, wall_temperature: f64) -> f64 {
    let film = 0.5 * (fluid_temperature + wall_temperature);
    let rho = 1000.0 - 0.5 * (film - 300.0);
    let beta = 0.5 / rho;
    let nu = 0.002 / rho;
    let alpha = 0.5 / (rho * 2000.0);
    let rayleigh = (9.81 * beta * (wall_temperature - fluid_temperature) * 0.1_f64.powi(3)
        / (nu * alpha)).abs();
    let prandtl = 8.0_f64;
    let nusselt = (0.825 + 0.387 * rayleigh.powf(1.0 / 6.0)
        / (1.0 + (0.492 / prandtl).powf(9.0 / 16.0)).powf(8.0 / 27.0)).powi(2);
    0.002 * nusselt * 0.5 / 0.1
}

#[test]
fn g1_natural_reservoir_closes_native_quasi_steady_storage_at_the_actual_endpoint() {
    context(|cx| {
        let card = buoyant_card();
        let (law, caloric) = source_pair(&card);
        let graph = graph();
        let laws = laws(2300.0);
        for wall_temperature in [280.0, 400.0] {
            let initial = with_reservoir(&caloric, [300.0; 6], Some(wall(wall_temperature)));
            let inventory = bind_inventory(cx, &graph, &initial);
            let owner = LiquidGraphSources::try_new(
                cx, &law, &caloric, &initial,
                SourceGraphBinding {
                    graph: &graph, inventory: &inventory, laws: &laws,
                    // Gravity is transverse to every original pipe section.
                    acceleration_m_s2: [9.81, 0.0, 0.0],
                }, caps(),
            ).unwrap();
            let problem = QuasiSteadyGraphProblem {
                topology: &graph, node_cells: inventory.node_cells(),
                reference_density: Density::new(1000.0), maximum_connections: 32,
            };
            let controls = QuasiSteadyGraphControls {
                solver: NewtonKrylovConfig {
                    absolute_tolerance: 1e-9, relative_tolerance: 0.0,
                    ..NewtonKrylovConfig::default()
                },
                maximum_cells: 16, maximum_newton_iterations: 32,
                maximum_endpoint_evaluations: 2048, energy_scale_j: 1.0,
                finite_difference_step: 1e-4, cell_tolerance_j: 1e-6,
                total_tolerance_j: 1e-6, maximum_neglected_power_w: 100.0,
                maximum_neglected_energy_j: 100.0,
            };
            let bindings = [GraphConvectiveReservoir {
                fluid_cell: CellId(2), characteristic_length: Length::new(0.1),
                property_temperature: ConvectionPropertyTemperature::FilmMean,
                law: GraphConvectionLaw::Natural { geometry: NaturalConvectionGeometry::VerticalPlate },
            }];
            let dt = 0.25;
            let (next, report) = advance_quasi_steady_graph_distributed(
                cx, &problem, &inventory, &initial, dt, controls, &[0.0; 6],
                |cx, candidate| {
                    let resolved = owner.resolve(cx, candidate, &[0.0; 3]).map_err(|e| e.to_string())?;
                    let hydraulics = resolved.graph.solve(
                        &[FixedPressure { node: 0, pressure: Pressure::new(0.0) }],
                        hydraulic_controls(), cx,
                    ).map_err(|e| e.to_string())?;
                    let heat = owner.evaluate_quasi_steady_reservoirs(
                        cx, candidate, &resolved.sources, &hydraulics, dt, &bindings, 1,
                    ).map_err(|e| e.to_string())?;
                    let mut boundaries = BoundaryConductances::default();
                    heat.convection.apply_to(&mut boundaries, || cx.checkpoint().is_err())
                        .map_err(|e| e.to_string())?;
                    Ok::<_, String>(QuasiSteadyGraphResponse {
                        hydraulics, mass_boundaries: vec![],
                        heat: BoundaryStepCoefficients {
                            bulk: StepCoefficients {
                                sources_w: BTreeMap::new(), conductances_w_k: BTreeMap::new(),
                                payload: heat,
                            }, boundaries,
                        },
                    })
                },
            ).unwrap();
            // One half-metre section contains 0.5 kg; all other cells remain
            // at rest and at their original temperature. Solve its independent
            // nonlinear backward-Euler balance, retaining the actual film law.
            let mut lower = wall_temperature.min(300.0);
            let mut upper = wall_temperature.max(300.0);
            for _ in 0..80 {
                let t = 0.5 * (lower + upper);
                let residual = 1000.0 * (t - 300.0)
                    - dt * natural_conductance(t, wall_temperature) * (wall_temperature - t);
                if residual > 0.0 { upper = t; } else { lower = t; }
            }
            let reference = 0.5 * (lower + upper);
            let actual = next.temperature_k(CellId(2)).unwrap();
            assert!((actual - reference).abs() < 1e-8);
            let sample = &report.heat.convection.samples()[0];
            near(sample.fluid_temperature_k, actual);
            near(sample.property_temperature_k, 0.5 * (actual + wall_temperature));
            near(sample.conductance_w_k, natural_conductance(actual, wall_temperature));
            assert_eq!(sample.heat_into_fluid_w.is_sign_positive(), wall_temperature > 300.0);
            let stored = next.energy_j() - initial.energy_j();
            assert!((stored - dt * sample.heat_into_fluid_w).abs() < 2e-6);
            assert!((report.thermal.external_heat_j - stored).abs() < 2e-6);
            assert!(report.hydraulics.branches.iter().all(|b| b.flow.value().abs() < 1e-14));
            for cell in [CellId(0), CellId(1), CellId(3), CellId(4), CellId(5)] {
                assert!((next.temperature_k(cell).unwrap() - 300.0).abs() < 1e-9);
            }
            near(next.time_s(), dt);
            near(initial.time_s(), 0.0);
            near(initial.temperature_k(CellId(2)).unwrap(), 300.0);
        }
    });
}

#[test]
fn g4_reservoir_identity_ownership_and_crossing_flow_keep_original_refusals() {
    context(|cx| {
        let card = card(true, 1.0);
        let (law, caloric) = source_pair(&card);
        let initial = with_reservoir(&caloric, [300.0; 6], Some(wall(400.0)));
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(1e6);
        let owner = LiquidGraphSources::try_new(
            cx, &law, &caloric, &initial,
            SourceGraphBinding {
                graph: &graph, inventory: &inventory, laws: &laws,
                acceleration_m_s2: [0.0; 3],
            }, caps(),
        ).unwrap();
        let rates = [1e-5, 0.0, 1e-5];
        let resolved = owner.resolve(cx, &initial, &rates).unwrap();
        let report = native_report(cx, &resolved.graph, &inventory, rates, [1.0, 0.0, 1.0], rates, 0.25);
        owner.evaluate_reservoirs(cx, &initial, &resolved.sources, &report, &[forced()], 1).unwrap();
        for replacement in [None, Some(Reservoir { conductance_w_k: 1.0, ..wall(400.0) })] {
            let changed = with_reservoir(&caloric, [300.0; 6], replacement);
            assert!(owner.evaluate_reservoirs(
                cx, &changed, &resolved.sources, &report, &[forced()], 1,
            ).is_err());
        }
        assert!(owner.evaluate_reservoirs(
            cx, &initial, &resolved.sources, &report, &[forced(), forced()], 2,
        ).is_err());
        let stale = with_reservoir(&caloric, [300.0, 300.0, 301.0, 300.0, 300.0, 300.0], Some(wall(400.0)));
        assert!(matches!(owner.evaluate_reservoirs(
            cx, &stale, &resolved.sources, &report, &[forced()], 1,
        ), Err(GraphReservoirConvectionError::Source(SourceGraphError::Invalid(_)))));
        let unsupported_film = with_reservoir(&caloric, [300.0; 6], Some(wall(650.0)));
        assert!(matches!(owner.evaluate_reservoirs(
            cx, &unsupported_film, &resolved.sources, &report, &[forced()], 1,
        ), Err(GraphReservoirConvectionError::Convection(_))));

        // Original inertance is 1e6 Pa s2/m3 per full branch and drag is
        // 1e5 Pa s/m3. This trajectory crosses from -0.0006 to +0.0012.
        // Mean, old and new |Re| all admit Gnielinski, but crossing zero
        // does not. The hydraulic section law deliberately covers the range.
        let average = [0.0003, 0.0, 0.0003];
        let crossing_sources = owner.resolve(cx, &initial, &average).unwrap();
        let crossing = native_report(
            cx, &crossing_sources.graph, &inventory,
            [-0.0006, 0.0, -0.0006], [18030.0, 0.0, 18030.0], average, 0.1,
        );
        let turbulent = GraphConvectiveReservoir {
            property_temperature: ConvectionPropertyTemperature::Reservoir,
            law: GraphConvectionLaw::Forced {
                branch: 0, correlation: CorrelationId::Gnielinski,
                settings: developed_settings(),
            }, ..forced()
        };
        assert!(matches!(owner.evaluate_reservoirs(
            cx, &initial, &crossing_sources.sources, &crossing, &[turbulent], 1,
        ), Err(GraphReservoirConvectionError::Envelope { .. })));
    });
}
