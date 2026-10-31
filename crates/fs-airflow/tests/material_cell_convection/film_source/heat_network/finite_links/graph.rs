//! Original material equations and finite coolant share one thermal endpoint.
use super::*;
use fs_couple::constitutive_graph::enthalpy::{
    GraphHeatPort, GraphTemperature, GraphThermalError, GraphThermalInput, GraphThermalSolve,
    ThermalCoordinate, heat_assembly_candidate,
};
use fs_couple::constitutive_graph::{
    MaterialEquation, MaterialGraphSolve, MaterialUnknown, MaterialValue, PortKey,
};
use fs_material::conductor::{graph::MaterialOhmicNode, resolve_conductor_state};
use fs_material::graph::ConstitutiveGraph;
use fs_qty::Dims;
use std::collections::BTreeMap;

const VOLTAGE: f64 = 2.0;
const CURRENT_DIMS: Dims = Dims([0, 0, 0, 0, 1, 0]);
const VOLTAGE_DIMS: Dims = Dims([2, 1, -3, 0, -1, 0]);

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
}

// This is one fixed manufactured specimen; its lot is identity metadata,
// not an independently evolving physical coordinate of the graph law.
fn graph_material(resistance_scale: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = body_material();
    let mut claims = ClaimSet::new();
    for (_, source) in original.claims().claims_ordered() {
        let mut claim = source.clone();
        claim.validity = ValidityDomain::unconstrained().with_quantity(
            "T",
            temperature_quantity(),
            280.0,
            320.0,
        );
        if claim.key.name() == "electrical_resistivity" {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                for (_, value) in knots {
                    *value *= resistance_scale;
                }
            }
        }
        claims.insert_claim(claim).unwrap();
    }
    let mut id = original.id().clone();
    id.process.push_str("; fixed graph specimen lot 7");
    let card = MaterialCard::assemble(id, claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn key(port: &str) -> PortKey {
    ("heater".into(), port.into())
}

fn material_graph(card: &MaterialCard) -> ConstitutiveGraph {
    let conductor = resolve_conductor_state(
        card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let node = MaterialOhmicNode::try_new(card, conductor, "T", [280.0, 320.0], 0.1, 1e-5).unwrap();
    let mut graph = ConstitutiveGraph::new();
    graph.add_node("heater", Box::new(node)).unwrap();
    graph
}

struct Feedback {
    fixed: BTreeMap<PortKey, f64>,
    unknowns: [MaterialUnknown; 1],
    equations: [MaterialEquation; 1],
    temperatures: [GraphTemperature; 1],
    heat: [GraphHeatPort; 1],
}
impl Feedback {
    fn new() -> Self {
        Self {
            fixed: BTreeMap::new(),
            unknowns: [MaterialUnknown {
                input: key("current"),
                dims: CURRENT_DIMS,
                scale_si: 1.0,
            }],
            equations: [MaterialEquation {
                terms: vec![(1.0, MaterialValue::Output(key("voltage")))],
                dims: VOLTAGE_DIMS,
                target_si: VOLTAGE,
                scale_si: VOLTAGE,
            }],
            temperatures: [GraphTemperature {
                input: key("temperature_change"),
                cell: BODY,
                coordinate: ThermalCoordinate::ChangeFrom(300.0),
            }],
            heat: [GraphHeatPort {
                output: key("heat"),
                cell: BODY,
            }],
        }
    }
    fn input<'a>(
        &'a self,
        graph: &'a ConstitutiveGraph,
        state: &'a [f64],
    ) -> GraphThermalInput<'a> {
        GraphThermalInput {
            graph,
            version: graph.state_schema().version(),
            state,
            fixed: &self.fixed,
            unknowns: &self.unknowns,
            equations: &self.equations,
            initial_coordinates: &[0.0],
            temperatures: &self.temperatures,
            durations: &[],
            heat: &self.heat,
        }
    }
}

fn solve_controls() -> GraphThermalSolve {
    GraphThermalSolve {
        maximum_cells: 2,
        maximum_thermal_trials: 400,
        maximum_graph_evaluations: 4000,
        maximum_iterations: 24,
        difference_step_j: 1e-4,
        cell_tolerance_j: 1e-8,
        total_tolerance_j: 1e-8,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 2,
            ..NewtonKrylovConfig::default()
        },
        graph: MaterialGraphSolve {
            maximum_unknowns: 1,
            maximum_terms: 8,
            maximum_evaluations: 100,
            maximum_iterations: 16,
            equation_tolerance: 1e-11,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-12,
                relative_tolerance: 1e-12,
                linear_restart: 1,
                ..NewtonKrylovConfig::default()
            },
        },
    }
}

// Independent voltage constraint, total heat and flat-plate equations.
// Eliminate coolant T by conservation, then solve the remaining wall balance.
fn reference(before: [f64; 2], viscosity: f64, cp: f64, resistance_scale: f64) -> [f64; 5] {
    let [wall0, fluid0] = before;
    let evaluate = |wall: f64| {
        let resistance = resistance_scale * (1.0 + 0.02 * (wall - 300.0));
        let current = VOLTAGE / resistance;
        let heat = DT * VOLTAGE * current;
        let fluid = fluid0 + (heat - MASS * CP * (wall - wall0)) / (MASS * cp);
        let film = 0.5 * (wall + fluid);
        let k = 0.5 + 0.002 * (film - 300.0);
        let mu = viscosity * 0.001 * (1.0 - 0.01 * (film - 300.0));
        let h = 0.664 * (1000.0 * 0.1 / mu).sqrt() * (mu * cp / k).cbrt() * k / 0.1;
        let power = SURFACE * h * (wall - fluid);
        (
            MASS * CP * (wall - wall0) - heat + DT * power,
            [wall, fluid, current, heat, power],
        )
    };
    let (mut low, mut high) = (280.0, 320.0);
    assert!(evaluate(low).0 < 0.0 && evaluate(high).0 > 0.0);
    for _ in 0..100 {
        let wall = 0.5 * (low + high);
        if evaluate(wall).0 > 0.0 {
            high = wall;
        } else {
            low = wall;
        }
    }
    let (residual, result) = evaluate(0.5 * (low + high));
    assert!(residual.abs() < 1e-10);
    result
}

#[test]
fn g1_g3_material_graph_current_and_finite_coolant_close_the_same_feedback() {
    let run = |viscosity, cp, resistance_scale| {
        let (card, curve) = graph_material(resistance_scale);
        let (fluid_card, fluid_curve) = fluid_material(viscosity, cp);
        let context = point(310.0);
        let fluid_context = fluid_point(290.0);
        let thermal = initial_pair(&curve, &fluid_curve, cp, false, 0.0);
        let graph = material_graph(&card);
        let feedback = Feedback::new();
        let fluid = liquid_law(&fluid_card);
        let cells = [
            MaterialCell {
                cell: BODY,
                material_id: MaterialId(0),
                curve: curve.curve(),
                point: ThermalMaterialPoint {
                    material: &card,
                    context: &context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: None,
            },
            MaterialCell {
                cell: FLUID,
                material_id: MaterialId(1),
                curve: fluid_curve.curve(),
                point: ThermalMaterialPoint {
                    material: &fluid_card,
                    context: &fluid_context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: None,
            },
        ];
        // The graph's original Ohmic node owns electrical work and Joule heat.
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &[],
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
        };
        let heat = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let links = [binding()];
        let initial = evaluate_source_link_convection(
            &thermal,
            &links,
            1,
            |_, _, t| fluid.at_temperature(t),
            || false,
        )
        .unwrap();
        let saved_cells = thermal.cells().clone();
        let saved_links = thermal.links().clone();
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let mut result = [0.0; 2];
        pool.scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 191,
                    kernel_id: 41,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            assert!(heat.prepare(&cx, &thermal).is_err());
            let prepared = heat
                .prepare_with_link_conductances(&cx, &thermal, initial.conductances_w_k())
                .unwrap();
            let network = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    fluid.at_temperature(t)
                },
            };
            let assembly = network.with_links(&links, 1, |_, _, t| fluid.at_temperature(t));
            let solve = solve_controls();
            let mut accepted = thermal.clone();
            let mut graph_state = Vec::new();
            let mut state = [310.0, 290.0];
            for step in 1..=2 {
                let before = accepted.cells().clone();
                let history = graph_state.clone();
                let expected = reference(state, viscosity, cp, resistance_scale);
                let next = heat_assembly_candidate(
                    &assembly,
                    &cx,
                    &accepted,
                    f64::from(step) * DT,
                    feedback.input(&graph, &graph_state),
                    solve,
                    &[0.0; 2],
                )
                .unwrap();
                let report = next.report();
                let current = next.graph().external[&key("current")];
                near(
                    next.thermal().temperature_k(BODY).unwrap(),
                    expected[0],
                    2e-8,
                );
                near(
                    next.thermal().temperature_k(FLUID).unwrap(),
                    expected[1],
                    2e-8,
                );
                near(current, expected[2], 1e-9);
                near(
                    next.graph().external[&key("temperature_change")],
                    expected[0] - 300.0,
                    2e-8,
                );
                near(next.graph().graph.outputs[&key("voltage")], VOLTAGE, 1e-10);
                near(
                    next.graph().graph.outputs[&key("heat")],
                    VOLTAGE * current,
                    1e-9,
                );
                near(report.graph_heat_j, expected[3], 1e-9);
                near(report.other_external_heat_j, 0.0, 1e-12);
                near(report.closed_radiation_net_j, 0.0, 1e-14);
                near(
                    next.thermal().energy_j() - accepted.energy_j(),
                    DT * VOLTAGE * current,
                    2e-9,
                );
                near(report.thermal.external_heat_j, report.graph_heat_j, 1e-9);
                assert!(report.maximum_cell_residual_j < solve.cell_tolerance_j);
                assert!(report.thermal_trials > 1 && report.graph_evaluations > 1);
                let finite = report.source.finite_convection.as_ref().unwrap();
                let sample = &finite.samples()[0];
                near(sample.heat_into_fluid_w, expected[4], 2e-7);
                near(MASS * cp * (expected[1] - state[1]), DT * expected[4], 1e-9);
                assert_eq!((sample.wall_cell, sample.fluid_cell), (BODY, FLUID));
                assert_eq!(
                    sample.material.material().card_identity(),
                    fluid_card.content_hash()
                );
                assert_eq!(sample.material.material().query().axes()["p"], 100000.0);
                assert_eq!(
                    sample.material.material().query().axes()["solute_mass_fraction"],
                    0.2
                );
                for cell in &report.source.heat.material.cells {
                    near(
                        cell.phase.temperature_k(),
                        next.thermal().temperature_k(cell.cell).unwrap(),
                        1e-12,
                    );
                }
                assert_eq!(accepted.cells(), &before);
                assert_eq!(graph_state, history);
                assert_eq!(next.thermal().links(), &saved_links);
                assert!(next.graph().graph.next_state.is_empty()); // Original Ohmic law is memoryless.
                result = [next.thermal().temperature_k(BODY).unwrap(), current];
                state = [expected[0], expected[1]];
                let (thermal_next, graph_next, _) = next.into_parts();
                accepted = thermal_next;
                graph_state = graph_next.graph.next_state;
            }
            if viscosity == 1.0 && cp == 4000.0 && resistance_scale == 1.0 {
                let mut insufficient = solve;
                insufficient.maximum_thermal_trials = 1;
                assert!(matches!(
                    heat_assembly_candidate(
                        &assembly,
                        &cx,
                        &thermal,
                        DT,
                        feedback.input(&graph, &[]),
                        insufficient,
                        &[0.0; 2]
                    ),
                    Err(GraphThermalError::ThermalBudget)
                ));
                let queries = std::cell::Cell::new(0);
                let refused = SourceConvectionHeatNetwork {
                    heat: &prepared,
                    cells: &[],
                    maximum_cells: 0,
                    resolve_fluid: |_: &SourceConvectiveCell,
                                    _: &DepositionSimulation<'_>,
                                    t: Temperature| {
                        Ok::<_, &'static str>(fluid.at_temperature(t).unwrap())
                    },
                };
                let failing = refused.with_links(&links, 1, |_, _, t| {
                    queries.set(queries.get() + 1);
                    if queries.get() > 1 {
                        Err("graph coolant trial refused")
                    } else {
                        Ok(fluid.at_temperature(t).unwrap())
                    }
                });
                assert!(matches!(
                    heat_assembly_candidate(
                        &failing,
                        &cx,
                        &thermal,
                        DT,
                        feedback.input(&graph, &[]),
                        solve,
                        &[0.0; 2]
                    ),
                    Err(GraphThermalError::Source(ConvectionHeatError::Convection(
                        SourceCellConvectionError::Source {
                            cell: FLUID,
                            source: "graph coolant trial refused",
                            ..
                        }
                    )))
                ));
                assert!(queries.get() > 1);
                let retry = heat_assembly_candidate(
                    &assembly,
                    &cx,
                    &thermal,
                    DT,
                    feedback.input(&graph, &[]),
                    solve,
                    &[0.0; 2],
                )
                .unwrap();
                let repeated = heat_assembly_candidate(
                    &assembly,
                    &cx,
                    &thermal,
                    DT,
                    feedback.input(&graph, &[]),
                    solve,
                    &[0.0; 2],
                )
                .unwrap();
                assert_eq!(retry.thermal().cells(), repeated.thermal().cells());
                assert_eq!(retry.graph().graph.outputs, repeated.graph().graph.outputs);
                gate.request();
                assert!(matches!(
                    heat_assembly_candidate(
                        &assembly,
                        &cx,
                        &thermal,
                        DT,
                        feedback.input(&graph, &[]),
                        solve,
                        &[0.0; 2]
                    ),
                    Err(GraphThermalError::Cancelled)
                ));
            }
            assert_eq!(thermal.cells(), &saved_cells);
            assert_eq!(thermal.links(), &saved_links);
            assert_eq!(thermal.time_s(), 0.0);
        });
        result
    };
    let baseline = run(1.0, 4000.0, 1.0);
    let viscous = run(4.0, 4000.0, 1.0);
    let higher_capacity = run(1.0, 8000.0, 1.0);
    let higher_resistivity = run(1.0, 4000.0, 2.0);
    assert!(viscous[0] > baseline[0] && viscous[1] < baseline[1]);
    assert!(higher_capacity[0] < baseline[0] && higher_capacity[1] > baseline[1]);
    assert!(higher_resistivity[0] < baseline[0] && higher_resistivity[1] < baseline[1]);
}
