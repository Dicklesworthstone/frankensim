//! The same sourced thermocouple through the public constitutive graph.
//! Native circuit fixtures remain the common physical description. A separate
//! four-cell linear heat equation supplies the prescribed-current G1 oracle.
use super::*;
use fs_couple::circuit_enthalpy::contacts::heat::{
    HeatNetworkError, HeatNetworkSample, MaterialHeatNetwork,
};
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
use fs_couple::constitutive_graph::enthalpy::*;
use fs_couple::constitutive_graph::*;
use fs_exec::{CancelGate, Cx};
use fs_material::conductor::{graph::MaterialOhmicNode, resolve_conductor_state};
use fs_material::graph::ConstitutiveGraph;
use fs_material::thermoelectric::graph::MaterialSeebeckNode;
use fs_qty::Dims;
use fs_solver::NewtonKrylovConfig;
use std::collections::BTreeMap;

const NODES: [&str; 4] = ["kelvin_a", "ohmic_a", "kelvin_b", "ohmic_b"];
const CURRENT: Dims = Dims([0, 0, 0, 0, 1, 0]);
const VOLTAGE: Dims = Dims([2, 1, -3, 0, -1, 0]);

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature,
        )
        .unwrap()
}

// These authored fixtures describe one fixed specimen. Its lot is specimen
// metadata, not an independently varying physical coordinate owned by a graph
// input. Both the graph and native circuit use these SAME cards and curves.
fn source_value(
    value: PropertyValue,
    interpolation: InterpolationPolicy,
    resistance: f64,
    conductivity: f64,
) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = super::source_value(value, interpolation, resistance, conductivity);
    let mut claims = fs_matdb::ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        let mut validity = ValidityDomain::unconstrained();
        for (axis, &(low, high)) in claim.validity.bounds() {
            if axis != "lot" {
                validity = if let Some(quantity) = claim.validity.axis_quantities().get(axis) {
                    validity.with_quantity(axis, *quantity, low, high)
                } else {
                    validity.with(axis, low, high)
                };
            }
        }
        claim.validity = validity;
        claim
            .provenance
            .source
            .push_str("; fixed authored specimen lot 7");
        claims.insert_claim(claim).unwrap();
    }
    let mut identity = original.id().clone();
    identity.process.push_str("; fixed specimen lot 7");
    let card = MaterialCard::assemble(identity, claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, caloric)
}

fn fixture(
    s1: f64,
    s2: f64,
    r2: f64,
    conductivity: f64,
) -> ([MaterialCard; 4], [ResolvedHeatCapacityCurve; 4]) {
    let [(a, ca), (b, cb), (c, cc), (d, cd)] = [
        (0.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
        (s1, 2.0, conductivity),
        (s2, r2, conductivity),
    ]
    .map(|(s, r, k)| {
        source_value(
            PropertyValue::Scalar {
                value: s,
                dims: SEEBECK_DIMS,
            },
            InterpolationPolicy::ConstantWithinValidity,
            r,
            k,
        )
    });
    ([a, b, c, d], [ca, cb, cc, cd])
}

fn port(node: &str, name: &str) -> PortKey {
    (node.into(), name.into())
}
fn with_cx<R>(run: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 61,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn material_graph(cards: &[MaterialCard; 4]) -> ConstitutiveGraph {
    let mut graph = ConstitutiveGraph::new();
    let context = point(300.0);
    for (kelvin, ohmic, cell) in [(NODES[0], NODES[1], 2), (NODES[2], NODES[3], 3)] {
        graph
            .add_node(
                kelvin,
                Box::new(
                    MaterialSeebeckNode::try_new(
                        &cards[cell],
                        &context,
                        "T",
                        [200.0, 2000.0],
                        MaterialPropertySelection::SingleClaimOnly,
                    )
                    .unwrap(),
                ),
            )
            .unwrap();
        let source = resolve_conductor_state(
            &cards[cell],
            &context,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        graph
            .add_node(
                ohmic,
                Box::new(
                    MaterialOhmicNode::try_new(
                        &cards[cell],
                        source,
                        "T",
                        [200.0, 2000.0],
                        1.0,
                        1.0,
                    )
                    .unwrap(),
                ),
            )
            .unwrap();
    }
    graph
}

fn controls() -> GraphThermalSolve {
    GraphThermalSolve {
        maximum_cells: 4,
        maximum_thermal_trials: 500,
        maximum_graph_evaluations: 10000,
        maximum_iterations: 24,
        difference_step_j: 1e-3,
        cell_tolerance_j: 1e-7,
        total_tolerance_j: 1e-7,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
        graph: MaterialGraphSolve {
            maximum_unknowns: 4,
            maximum_terms: 24,
            maximum_evaluations: 200,
            maximum_iterations: 16,
            equation_tolerance: 1e-10,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-12,
                relative_tolerance: 1e-12,
                linear_restart: 8,
                ..NewtonKrylovConfig::default()
            },
        },
    }
}

fn graph_endpoint<'m>(
    cx: &Cx<'_>,
    cards: &[MaterialCard; 4],
    curves: &[ResolvedHeatCapacityCurve; 4],
    initial: &DepositionSimulation<'m>,
    target: f64,
    current: Option<f64>,
    voltage: f64,
) -> Result<GraphThermalCandidate<'m, HeatNetworkSample>, GraphThermalError<HeatNetworkError>> {
    let graph = material_graph(cards);
    let mut fixed = BTreeMap::new();
    let mut unknowns = Vec::new();
    let mut equations = Vec::new();
    if let Some(current) = current {
        for node in NODES {
            fixed.insert(port(node, "current"), current);
        }
    } else {
        for node in NODES {
            unknowns.push(MaterialUnknown {
                input: port(node, "current"),
                dims: CURRENT,
                scale_si: 1.0,
            });
        }
        for node in &NODES[1..] {
            equations.push(MaterialEquation {
                terms: vec![
                    (1.0, MaterialValue::Input(port(NODES[0], "current"))),
                    (-1.0, MaterialValue::Input(port(node, "current"))),
                ],
                dims: CURRENT,
                target_si: 0.0,
                scale_si: 1.0,
            });
        }
        equations.push(MaterialEquation {
            terms: NODES
                .into_iter()
                .map(|node| (1.0, MaterialValue::Output(port(node, "voltage"))))
                .collect(),
            dims: VOLTAGE,
            target_si: voltage,
            scale_si: 1.0,
        });
    }
    let guess = vec![0.0; unknowns.len()];
    let mut temperatures = Vec::new();
    let mut heat = Vec::new();
    for (kelvin, ohmic, ids) in [
        (NODES[0], NODES[1], [0, 2, 1]),
        (NODES[2], NODES[3], [1, 3, 0]),
    ] {
        for (name, cell) in [
            "temperature_change_a",
            "temperature_change_bulk",
            "temperature_change_b",
        ]
        .into_iter()
        .zip(ids)
        {
            temperatures.push(GraphTemperature {
                input: port(kelvin, name),
                cell: CellId(cell),
                coordinate: ThermalCoordinate::ChangeFrom(300.0),
            });
        }
        for (name, cell) in ["heat_a", "heat_bulk", "heat_b"].into_iter().zip(ids) {
            heat.push(GraphHeatPort {
                output: port(kelvin, name),
                cell: CellId(cell),
            });
        }
        temperatures.push(GraphTemperature {
            input: port(ohmic, "temperature_change"),
            cell: CellId(ids[1]),
            coordinate: ThermalCoordinate::ChangeFrom(300.0),
        });
        heat.push(GraphHeatPort {
            output: port(ohmic, "heat"),
            cell: CellId(ids[1]),
        });
    }
    let context = point(300.0);
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(cards, curves, &context, &conductivity);
    let links = links();
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &links,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat_network = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    heat_network
        .prepare(cx, initial)
        .map_err(GraphThermalError::Source)?
        .graph_step_candidate(
            cx,
            initial,
            target,
            GraphThermalInput {
                graph: &graph,
                version: graph.state_schema().version(),
                state: &[],
                fixed: &fixed,
                unknowns: &unknowns,
                equations: &equations,
                initial_coordinates: &guess,
                temperatures: &temperatures,
                durations: &[],
                heat: &heat,
            },
            controls(),
            &[0.0; 4],
        )
}

fn native_endpoint<'m>(
    cards: &[MaterialCard; 4],
    curves: &[ResolvedHeatCapacityCurve; 4],
    initial: &DepositionSimulation<'m>,
    target: f64,
    current: Option<f64>,
    voltage: f64,
) -> thermoelectric::ThermoelectricCandidate<'m> {
    let context = point(300.0);
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(cards, curves, &context, &conductivity);
    let resistors = resistors();
    let links = links();
    let material = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let drive = current.unwrap_or(voltage);
    let electrical = circuit(&material, initial, current.is_some(), drive);
    thermoelectric::step_candidate(
        &material,
        initial,
        target,
        &electrical,
        &[0.0, 0.0, drive],
        &legs(),
        &BASES,
        config(),
        &[0.0; 4],
        || false,
    )
    .unwrap()
}

// Independent constant-property BE equations for the four actual physical cells.
// Q_cold=-I*(Sa-Sb)*T_cold, Q_hot=+I*(Sa-Sb)*T_hot, and each bulk
// receives its own I²R. Fourier transfer is an explicit symmetric graph Laplacian.
fn linear_reference(
    old: [f64; 4],
    dt: f64,
    current: f64,
    seebeck_difference: f64,
    resistance_b: f64,
    conductivity: f64,
) -> [f64; 4] {
    let mut a = [[0.0; 5]; 4];
    for i in 0..4 {
        a[i][i] = 10.0 / dt;
        a[i][4] = 10.0 / dt * old[i];
    }
    a[0][0] += current * seebeck_difference;
    a[1][1] -= current * seebeck_difference;
    a[2][4] += current * current * 2.0;
    a[3][4] += current * current * resistance_b;
    let g = 1.0 / (0.5 + 0.5 / conductivity);
    for [i, j] in [[0, 2], [2, 1], [1, 3], [3, 0]] {
        a[i][i] += g;
        a[j][j] += g;
        a[i][j] -= g;
        a[j][i] -= g;
    }
    for column in 0..4 {
        let pivot = a[column][column];
        assert!(pivot.abs() > 1.0);
        for j in column..5 {
            a[column][j] /= pivot;
        }
        for i in 0..4 {
            if i != column {
                let multiplier = a[i][column];
                for j in column..5 {
                    a[i][j] -= multiplier * a[column][j];
                }
            }
        }
    }
    a.map(|row| row[4])
}

fn compare_temperatures(a: &DepositionSimulation<'_>, b: &DepositionSimulation<'_>) {
    for cell in 0..4 {
        near(
            a.temperature_k(CellId(cell)).unwrap(),
            b.temperature_k(CellId(cell)).unwrap(),
        );
    }
}
fn audit_transfer(candidate: &GraphThermalCandidate<'_, HeatNetworkSample>) {
    let outputs = &candidate.graph().graph.outputs;
    for (node, ids) in [(NODES[0], [0, 2, 1]), (NODES[2], [1, 3, 0])] {
        let heat = ["heat_a", "heat_bulk", "heat_b"].map(|p| outputs[&port(node, p)]);
        near(
            heat.iter().sum::<f64>() + outputs[&port(node, "electrical_supply")],
            0.0,
        );
        let entropy: f64 = heat
            .into_iter()
            .zip(ids)
            .map(|(q, id)| q / candidate.thermal().temperature_k(CellId(id)).unwrap())
            .sum();
        assert!(entropy.abs() < 1e-12, "Kelvin entropy transfer: {entropy}");
    }
    let joule = outputs[&port(NODES[1], "heat")] + outputs[&port(NODES[3], "heat")];
    near(candidate.graph().graph.total_dissipation, joule);
    assert!(joule >= 0.0);
}

#[test]
fn prescribed_current_cools_finite_junctions_and_matches_independent_heat_equations() {
    with_cx(|cx| {
        let (cards, curves) = fixture(0.02, -0.01, 3.0, 1.0);
        let initial = initial(&curves, [300.0; 4]);
        let original = initial.cells().clone();
        for current in [-0.5, 0.5] {
            let candidate =
                graph_endpoint(cx, &cards, &curves, &initial, 0.2, Some(current), 0.0).unwrap();
            let reference = linear_reference([300.0; 4], 0.2, current, 0.03, 3.0, 1.0);
            for (i, expected) in reference.into_iter().enumerate() {
                near(
                    candidate.thermal().temperature_k(CellId(i as u64)).unwrap(),
                    expected,
                );
            }
            let native = native_endpoint(&cards, &curves, &initial, 0.2, Some(current), 0.0);
            compare_temperatures(candidate.thermal(), native.thermal());
            near(
                candidate.report().graph_heat_j,
                native.report().external_electrical_work_j,
            );
            assert!(candidate.graph().coordinates.is_empty());
            assert_eq!(candidate.graph().nonlinear.iterations, 0);
            assert_eq!(candidate.graph().evaluations, 2);
            let cold = candidate.thermal().temperature_k(CellId(0)).unwrap();
            let hot = candidate.thermal().temperature_k(CellId(1)).unwrap();
            assert!((hot - cold) * current > 0.0);
            audit_transfer(&candidate);
        }
        assert_eq!(initial.cells(), &original);
        assert_eq!(initial.time_s(), 0.0);
    });
}

#[test]
fn material_swaps_close_generator_current_and_internal_kelvin_work_through_both_native_paths() {
    with_cx(|cx| {
        let mut answers = Vec::new();
        for (seebeck, resistance, conductivity) in [
            (0.02, 3.0, 1.0),
            (0.04, 3.0, 1.0),
            (0.02, 6.0, 1.0),
            (0.02, 3.0, 3.0),
        ] {
            let (cards, curves) = fixture(seebeck, -0.01, resistance, conductivity);
            let initial = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
            let candidate = graph_endpoint(cx, &cards, &curves, &initial, 0.1, None, 0.0).unwrap();
            let native = native_endpoint(&cards, &curves, &initial, 0.1, None, 0.0);
            compare_temperatures(candidate.thermal(), native.thermal());
            let current = candidate.graph().external[&port(NODES[0], "current")];
            for node in NODES {
                near(candidate.graph().external[&port(node, "current")], current);
            }
            near(
                current,
                native.report().coupled.source.legs[0].transfer.current_a,
            );
            let difference = candidate.thermal().temperature_k(CellId(0)).unwrap()
                - candidate.thermal().temperature_k(CellId(1)).unwrap();
            near(current, (seebeck + 0.01) * difference / (2.0 + resistance));
            near(candidate.report().graph_heat_j, 0.0);
            near(candidate.thermal().energy_j(), initial.energy_j());
            near(candidate.report().other_external_heat_j, 0.0);
            audit_transfer(&candidate);
            answers.push((current, difference));
        }
        assert!(answers[1].0 > answers[0].0);
        assert!(answers[2].0 < answers[0].0);
        assert!(answers[3].1 < answers[0].1);
    });
}

#[test]
fn variable_thermopower_drives_original_bulk_heat_on_continued_thermal_history() {
    with_cx(|cx| {
        let (mut cards, mut curves) = fixture(0.02, -0.01, 3.0, 1.0);
        let slope = 2e-5;
        (cards[2], curves[2]) = source_value(
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (200.0, 0.02 - 150.0 * slope),
                    (2000.0, 0.02 + 1650.0 * slope),
                ],
                dims: SEEBECK_DIMS,
            },
            InterpolationPolicy::LinearInside,
            2.0,
            1.0,
        );
        let initial = initial(&curves, [400.0, 300.0, 350.0, 350.0]);
        let first = graph_endpoint(cx, &cards, &curves, &initial, 0.15, Some(0.5), 0.0).unwrap();
        let second =
            graph_endpoint(cx, &cards, &curves, first.thermal(), 0.3, Some(0.5), 0.0).unwrap();
        let native = native_endpoint(&cards, &curves, first.thermal(), 0.3, Some(0.5), 0.0);
        compare_temperatures(second.thermal(), native.thermal());
        let [ta, tm, tb] = [0, 2, 1].map(|id| second.thermal().temperature_k(CellId(id)).unwrap());
        let expected = 0.25 * slope * (ta - tb) * tm;
        near(
            second.graph().graph.outputs[&port(NODES[0], "heat_bulk")],
            expected,
        );
        assert!(expected > 0.0);
        audit_transfer(&second);
        near(
            second.thermal().energy_j() - first.thermal().energy_j(),
            second.report().graph_heat_j,
        );
        near(
            second.report().graph_heat_j,
            native.report().external_electrical_work_j,
        );
        assert_eq!(first.thermal().time_s(), 0.15);
        assert_eq!(second.thermal().time_s(), 0.3);
        assert_eq!(initial.time_s(), 0.0);
        let old = second.thermal().cells().clone();
        assert!(
            graph_endpoint(
                cx,
                &cards,
                &curves,
                second.thermal(),
                10.3,
                Some(100.0),
                0.0
            )
            .is_err()
        );
        assert_eq!(second.thermal().cells(), &old);
        assert_eq!(second.thermal().time_s(), 0.3);
    });
}
