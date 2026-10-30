//! Real material transport and native finite thermal storage, including latent heat.
#[path = "enthalpy/expansion.rs"]
mod expansion;
#[path = "enthalpy/spatial.rs"]
mod spatial;

use super::*;
use fs_couple::circuit_enthalpy::contacts::heat::MaterialHeatNetwork;
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
use fs_couple::circuit_enthalpy::material::{ConductionLink, MaterialCell, MaterialNetwork};
use fs_couple::constitutive_graph::enthalpy::*;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::conductor::{
    ELECTRICAL_RESISTIVITY_DIMS, graph::MaterialOhmicNode, resolve_conductor_state,
};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacitySourceQuery,
    ResolvedHeatCapacityCurve, SolidLiquidPhase,
};
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS,
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Density, Temperature};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use std::cell::RefCell;
use std::rc::Rc;

const HEATER: CellId = CellId(30);
const PASSIVE: CellId = CellId(10);
const RESISTANCE_SLOPE: f64 = 0.01;
const VOLTAGE: f64 = 10.0;

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new().with("T", temperature).unwrap()
}

// E, R and k remain independent DATA on one original card. Manufactured values
// make the endpoint equation analytical; these fixtures claim no measured alloy.
fn material(r0: f64, cp: f64, k_scale: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, low, high) in [
        (
            "electrical_resistivity",
            ELECTRICAL_RESISTIVITY_DIMS,
            r0 - 100.0 * RESISTANCE_SLOPE,
            r0 + 500.0 * RESISTANCE_SLOPE,
        ),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            cp,
            cp,
        ),
        (
            "thermal_conductivity",
            THERMAL_CONDUCTIVITY_DIMS,
            0.8 * k_scale,
            2.0 * k_scale,
        ),
        ("density", Density::DIMS, 1.0, 1.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, low), (800.0, high)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 200.0, 800.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured graph/enthalpy source data".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured".into(),
            phase: "solid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(800.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-7,
    )
    .unwrap();
    (card, curve)
}

fn ohmic_graph(card: &MaterialCard) -> ConstitutiveGraph {
    let original = resolve_conductor_state(
        card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let node = MaterialOhmicNode::try_new(card, original, "T", [200.0, 800.0], 1.0, 1.0).unwrap();
    let mut graph = ConstitutiveGraph::new();
    graph.add_node("resistor", Box::new(node)).unwrap();
    graph
}

fn boundary() -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: 0.0,
    }
}
fn birth(id: CellId, material: MaterialId, h: f64) -> TimedEvent {
    TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id,
            segment_id: id.0,
            body_id: id.0,
            material_id: material,
            mass_kg: 1.0,
            specific_h_j_kg: h,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [1.0, 0.0, 0.0],
                volume_m3: 1.0,
                surface_area_m2: 6.0,
            },
            boundary: boundary(),
        }),
    }
}
fn timeline<'a>(
    laws: &[(MaterialId, &'a dyn ThermalLaw)],
    events: Vec<TimedEvent>,
) -> DepositionSimulation<'a> {
    DepositionSimulation::new(
        events,
        laws,
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 4,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn thermal_controls() -> GraphThermalSolve {
    GraphThermalSolve {
        maximum_cells: 2,
        maximum_thermal_trials: 400,
        maximum_graph_evaluations: 4000,
        maximum_iterations: 24,
        difference_step_j: 1e-4,
        cell_tolerance_j: 1e-7,
        total_tolerance_j: 1e-7,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 4,
            ..NewtonKrylovConfig::default()
        },
        graph: controls(),
    }
}

struct Feedback {
    fixed: BTreeMap<PortKey, f64>,
    unknowns: Vec<MaterialUnknown>,
    equations: Vec<MaterialEquation>,
    temperatures: Vec<GraphTemperature>,
    durations: Vec<PortKey>,
    heat: Vec<GraphHeatPort>,
}
impl Feedback {
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
            durations: &self.durations,
            heat: &self.heat,
        }
    }
}
fn voltage_feedback() -> Feedback {
    Feedback {
        fixed: BTreeMap::new(),
        unknowns: vec![MaterialUnknown {
            input: key("resistor", "current"),
            dims: A,
            scale_si: 5.0,
        }],
        equations: vec![MaterialEquation {
            terms: vec![(1.0, MaterialValue::Output(key("resistor", "voltage")))],
            dims: V,
            target_si: VOLTAGE,
            scale_si: VOLTAGE,
        }],
        temperatures: vec![GraphTemperature {
            input: key("resistor", "temperature_change"),
            cell: HEATER,
            coordinate: ThermalCoordinate::ChangeFrom(300.0),
        }],
        durations: vec![],
        heat: vec![GraphHeatPort {
            output: key("resistor", "heat"),
            cell: HEATER,
        }],
    }
}
fn no_additional_heat(_: &DepositionSimulation<'_>) -> Result<StepCoefficients<()>, &'static str> {
    Ok(StepCoefficients {
        sources_w: BTreeMap::new(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    })
}

#[test]
fn sourced_resistance_and_heat_capacity_close_the_finite_cell_quadratic() {
    with_cx(|cx, _| {
        let mut endpoints = Vec::new();
        for (r0, cp) in [(2.0, 2.0), (4.0, 2.0), (2.0, 4.0)] {
            let (card, curve) = material(r0, cp, 1.0);
            assert!(!curve.receipts().is_empty());
            let initial = timeline(
                &[(MaterialId(0), curve.curve() as &dyn ThermalLaw)],
                vec![birth(HEATER, MaterialId(0), 100.0 * cp)],
            );
            let original_cells = initial.cells().clone();
            let graph = ohmic_graph(&card);
            let feedback = voltage_feedback();
            let dt = 0.5;
            let candidate = step_candidate(
                cx,
                &initial,
                dt,
                feedback.input(&graph, &[]),
                thermal_controls(),
                &[0.0],
                no_additional_heat,
            )
            .unwrap();
            // E*(R0 + a*E/C) = dt*V². This stable positive quadratic root
            // uses source constants, without replaying either native solver.
            let supply = dt * VOLTAGE * VOLTAGE;
            let expected_energy =
                2.0 * supply / (r0 + (r0 * r0 + 4.0 * RESISTANCE_SLOPE * supply / cp).sqrt());
            let temperature = candidate.thermal().temperature_k(HEATER).unwrap();
            let current = candidate.graph().external[&key("resistor", "current")];
            close(
                candidate.thermal().energy_j() - initial.energy_j(),
                expected_energy,
            );
            close(temperature, 300.0 + expected_energy / cp);
            close(
                current,
                VOLTAGE / (r0 + RESISTANCE_SLOPE * (temperature - 300.0)),
            );
            close(
                candidate.graph().external[&key("resistor", "temperature_change")],
                temperature - 300.0,
            );
            close(
                candidate.graph().graph.outputs[&key("resistor", "heat")],
                VOLTAGE * current,
            );
            close(candidate.report().graph_heat_j, expected_energy);
            close(candidate.report().thermal.external_heat_j, expected_energy);
            close(candidate.report().other_external_heat_j, 0.0);
            assert!(candidate.report().thermal_trials > 1);
            assert!(candidate.graph().graph.next_state.is_empty());
            assert_eq!(initial.cells(), &original_cells);
            assert_eq!(initial.time_s(), 0.0);
            endpoints.push((temperature, current));
        }
        assert!(endpoints[1].0 < endpoints[0].0 && endpoints[1].1 < endpoints[0].1);
        assert!(endpoints[2].0 < endpoints[0].0 && endpoints[2].1 > endpoints[0].1);
    });
}

#[test]
fn conductivity_data_changes_both_original_bodies_and_the_converged_current() {
    with_cx(|cx, _| {
        let run = |k_scale: f64| {
            let a = material(2.0, 2.0, k_scale);
            let b = material(3.0, 3.0, 0.7 * k_scale);
            let cards = [a.0, b.0];
            let curves = [a.1, b.1];
            let context = point(300.0);
            let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
            let ids = [HEATER, PASSIVE];
            let bindings = [0, 1].map(|i| MaterialCell {
                cell: ids[i],
                material_id: MaterialId(i as u64),
                curve: curves[i].curve(),
                point: ThermalMaterialPoint {
                    material: &cards[i],
                    context: &context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: Some((&conductivity, MaterialPropertySelection::SingleClaimOnly)),
            });
            let links = [ConductionLink {
                link: 7,
                distance_a_m: 0.5,
                distance_b_m: 0.5,
            }];
            let bulk = MaterialNetwork {
                cells: &bindings,
                resistors: &[],
                links: &links,
            };
            let electrical = ElectricalContactNetwork {
                bulk: &bulk,
                contacts: &[],
                electrical_bases: &[],
            };
            let network = MaterialHeatNetwork {
                electrical: &electrical,
                thermal_contacts: &[],
                radiation: None,
            };
            let initial = timeline(
                &[
                    (MaterialId(0), curves[0].curve() as &dyn ThermalLaw),
                    (MaterialId(1), curves[1].curve() as &dyn ThermalLaw),
                ],
                vec![
                    birth(HEATER, MaterialId(0), 2.0 * 140.0),
                    birth(PASSIVE, MaterialId(1), 3.0 * 100.0),
                    TimedEvent {
                        time_s: 0.0,
                        event: Event::Contact(ThermalLink {
                            id: 7,
                            a: HEATER,
                            b: PASSIVE,
                            area_m2: 1.0,
                            conductance_w_k: 999.0,
                        }),
                    },
                ],
            );
            let prepared = network.prepare(cx, &initial).unwrap();
            let graph = ohmic_graph(&cards[0]);
            let feedback = voltage_feedback();
            let dt = 0.25;
            let candidate = prepared
                .graph_step_candidate(
                    cx,
                    &initial,
                    dt,
                    feedback.input(&graph, &[]),
                    thermal_controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let temperatures = ids.map(|id| candidate.thermal().temperature_k(id).unwrap());
            let current = candidate.graph().external[&key("resistor", "current")];
            let report = candidate.report();
            close(
                candidate.thermal().energy_j() - initial.energy_j(),
                dt * VOLTAGE * current,
            );
            close(report.graph_heat_j, dt * VOLTAGE * current);
            close(report.thermal.external_heat_j, report.graph_heat_j);
            close(report.other_external_heat_j, 0.0);
            for sample in &report.source.material.cells {
                let i = if sample.cell == HEATER { 0 } else { 1 };
                close(sample.phase.temperature_k(), temperatures[i]);
                let source = sample.conductivity.as_ref().unwrap();
                close(source.query().axes()["T"], temperatures[i]);
                close(
                    source.property("thermal_conductivity").unwrap().value_si(),
                    k_scale
                        * if i == 0 { 1.0 } else { 0.7 }
                        * (1.0 + 0.002 * (temperatures[i] - 300.0)),
                );
            }
            let k = [
                k_scale * (1.0 + 0.002 * (temperatures[0] - 300.0)),
                0.7 * k_scale * (1.0 + 0.002 * (temperatures[1] - 300.0)),
            ];
            let conductance = 1.0 / (0.5 / k[0] + 0.5 / k[1]);
            close(report.source.material.links[0].conductance_w_k, conductance);
            close(
                3.0 * (temperatures[1] - 300.0),
                dt * conductance * (temperatures[0] - temperatures[1]),
            );
            close(
                current,
                VOLTAGE / (2.0 + RESISTANCE_SLOPE * (temperatures[0] - 300.0)),
            );
            close(
                candidate.graph().external[&key("resistor", "temperature_change")],
                temperatures[0] - 300.0,
            );
            assert_eq!(initial.time_s(), 0.0);
            (temperatures, current)
        };
        let low = run(0.2);
        let high = run(2.0);
        assert!(high.0[0] < low.0[0]);
        assert!(high.0[1] > low.0[1]);
        assert!(high.1 > low.1);
    });
}

struct HeatDose {
    declaration: NodeDeclaration,
    observed_history: Rc<RefCell<Vec<f64>>>,
}
impl LawNode for HeatDose {
    fn declaration(&self) -> &NodeDeclaration {
        &self.declaration
    }
    fn evaluate(&self, state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        // Test instrumentation only: physical memory is returned, never mutated.
        self.observed_history.borrow_mut().push(state[0]);
        Ok(NodeOutput {
            outputs: vec![inputs[0]],
            next_state: vec![state[0] + inputs[0] * inputs[1]],
            dissipation_rate: None,
        })
    }
    fn tangent(&self, _state: &[f64], _inputs: &[f64]) -> Option<Vec<f64>> {
        Some(vec![1.0, 0.0])
    }
}
fn dose_graph() -> (ConstitutiveGraph, Rc<RefCell<Vec<f64>>>) {
    let observed = Rc::new(RefCell::new(Vec::new()));
    let node = HeatDose {
        declaration: declaration(
            "manufactured-heat-dose",
            &[("power", W), ("dt", Dims([0, 0, 1, 0, 0, 0]))],
            &[("heat", W)],
            &["dose_j"],
        ),
        observed_history: observed.clone(),
    };
    let mut graph = ConstitutiveGraph::new();
    graph.add_node("dose", Box::new(node)).unwrap();
    observed.borrow_mut().clear();
    (graph, observed)
}
fn dose_feedback() -> Feedback {
    Feedback {
        fixed: BTreeMap::new(),
        unknowns: vec![MaterialUnknown {
            input: key("dose", "power"),
            dims: W,
            scale_si: 300.0,
        }],
        equations: vec![MaterialEquation {
            terms: vec![(1.0, MaterialValue::Output(key("dose", "heat")))],
            dims: W,
            target_si: 300.0,
            scale_si: 300.0,
        }],
        temperatures: vec![],
        durations: vec![key("dose", "dt")],
        heat: vec![GraphHeatPort {
            output: key("dose", "heat"),
            cell: HEATER,
        }],
    }
}
fn latent_curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        fs_blake3::hash_bytes(b"manufactured latent heat chart"),
        [
            (-100.0, 300.0, 0.0),
            (0.0, 400.0, 0.0),
            (100.0, 400.0, 1.0),
            (200.0, 500.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1.0,
        })
        .collect(),
    )
    .unwrap()
}

#[test]
fn latent_plateau_keeps_native_enthalpy_and_commits_graph_memory_once_per_interval() {
    with_cx(|cx, _| {
        let curve = latent_curve();
        let initial = timeline(
            &[(MaterialId(0), &curve as &dyn ThermalLaw)],
            vec![birth(HEATER, MaterialId(0), -20.0)],
        );
        let (graph, observed) = dose_graph();
        let feedback = dose_feedback();
        let original_history = [5.0];
        let first = step_candidate(
            cx,
            &initial,
            0.1,
            feedback.input(&graph, &original_history),
            thermal_controls(),
            &[0.0],
            no_additional_heat,
        )
        .unwrap();
        close(first.thermal().cells()[&HEATER].specific_h_j_kg, 10.0);
        close(first.thermal().temperature_k(HEATER).unwrap(), 400.0);
        close(
            curve
                .state_at_specific_enthalpy(first.thermal().cells()[&HEATER].specific_h_j_kg)
                .unwrap()
                .liquid_mass_fraction(),
            0.1,
        );
        close(first.graph().graph.next_state[0], 35.0);
        assert!(observed.borrow().len() > 1);
        assert!(observed.borrow().iter().all(|&history| history == 5.0));
        assert_eq!(original_history, [5.0]);
        assert_eq!(initial.cells()[&HEATER].specific_h_j_kg, -20.0);
        let (thermal, accepted_graph, first_report) = first.into_parts();
        let committed_history = accepted_graph.graph.next_state.clone();
        let committed_cells = thermal.cells().clone();
        observed.borrow_mut().clear();
        // The second interval is longer: dt must come from the accepted thermal
        // clock, not a stale graph input retained from the first step.
        let second = step_candidate(
            cx,
            &thermal,
            0.25,
            feedback.input(&graph, &accepted_graph.graph.next_state),
            thermal_controls(),
            &[0.0],
            no_additional_heat,
        )
        .unwrap();
        close(second.thermal().cells()[&HEATER].specific_h_j_kg, 55.0);
        close(second.thermal().temperature_k(HEATER).unwrap(), 400.0);
        close(
            curve
                .state_at_specific_enthalpy(second.thermal().cells()[&HEATER].specific_h_j_kg)
                .unwrap()
                .liquid_mass_fraction(),
            0.55,
        );
        close(second.graph().graph.next_state[0], 80.0);
        close(second.graph().external[&key("dose", "dt")], 0.15);
        assert!(observed.borrow().len() > 1);
        assert!(
            observed
                .borrow()
                .iter()
                .all(|&history| history == committed_history[0])
        );
        close(
            second.thermal().energy_j() - initial.energy_j(),
            first_report.graph_heat_j + second.report().graph_heat_j,
        );
        close(second.report().thermal.residual_j, 0.0);
        assert_eq!(thermal.cells(), &committed_cells);
        assert_eq!(accepted_graph.graph.next_state, committed_history);
    });
}

#[test]
fn source_refusal_duplicate_route_events_work_bounds_and_cancellation_keep_both_histories() {
    with_cx(|cx, gate| {
        let curve = latent_curve();
        let initial = timeline(
            &[(MaterialId(0), &curve as &dyn ThermalLaw)],
            vec![birth(HEATER, MaterialId(0), -20.0)],
        );
        let saved_cells = initial.cells().clone();
        let (graph, _) = dose_graph();
        let history = [5.0];
        let feedback = dose_feedback();
        let source_failure = step_candidate(
            cx,
            &initial,
            0.1,
            feedback.input(&graph, &history),
            thermal_controls(),
            &[0.0],
            |_| Err::<StepCoefficients<()>, _>("original heat-source refusal"),
        );
        assert!(matches!(
            source_failure,
            Err(GraphThermalError::Source("original heat-source refusal"))
        ));
        let mut duplicate = dose_feedback();
        duplicate.heat.push(duplicate.heat[0].clone());
        assert!(matches!(
            step_candidate(
                cx,
                &initial,
                0.1,
                duplicate.input(&graph, &history),
                thermal_controls(),
                &[0.0],
                no_additional_heat
            ),
            Err(GraphThermalError::Invalid(_))
        ));
        let crossing = timeline(
            &[(MaterialId(0), &curve as &dyn ThermalLaw)],
            vec![
                birth(HEATER, MaterialId(0), -20.0),
                TimedEvent {
                    time_s: 0.05,
                    event: Event::Boundary(
                        HEATER,
                        Boundary {
                            source_w: 10.0,
                            ..boundary()
                        },
                    ),
                },
            ],
        );
        assert!(matches!(
            step_candidate(
                cx,
                &crossing,
                0.1,
                feedback.input(&graph, &history),
                thermal_controls(),
                &[0.0],
                no_additional_heat
            ),
            Err(GraphThermalError::Thermal(DepositionError::InvalidInput(_)))
        ));
        let mut bounded = thermal_controls();
        bounded.maximum_thermal_trials = 1;
        assert!(matches!(
            step_candidate(
                cx,
                &initial,
                0.1,
                feedback.input(&graph, &history),
                bounded,
                &[0.0],
                no_additional_heat
            ),
            Err(GraphThermalError::ThermalBudget)
        ));
        bounded = thermal_controls();
        bounded.maximum_graph_evaluations = 1;
        assert!(matches!(
            step_candidate(
                cx,
                &initial,
                0.1,
                feedback.input(&graph, &history),
                bounded,
                &[0.0],
                no_additional_heat
            ),
            Err(GraphThermalError::Graph(
                MaterialGraphError::EvaluationBudget
            ))
        ));
        gate.request();
        assert!(matches!(
            step_candidate(
                cx,
                &initial,
                0.1,
                feedback.input(&graph, &history),
                thermal_controls(),
                &[0.0],
                no_additional_heat
            ),
            Err(GraphThermalError::Cancelled)
        ));
        assert_eq!(initial.cells(), &saved_cells);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(crossing.cells(), &saved_cells);
        assert_eq!(crossing.time_s(), 0.0);
        assert_eq!(history, [5.0]);
    });
}
