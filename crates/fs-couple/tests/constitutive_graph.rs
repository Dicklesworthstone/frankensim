//! Native material cards and closed electrothermal feedback through public residuals.
#[path = "constitutive_graph/enthalpy.rs"]
mod enthalpy;

use fs_couple::constitutive_graph::*;
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_material::graph::*;
use fs_material::hyper::cards::register_hyperelastic_laws;
use fs_qty::Dims;
use fs_solver::{NewtonKrylovConfig, NonlinearProblem};
use std::collections::BTreeMap;

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const K: Dims = Dims([0, 0, 0, 1, 0, 0]);
const A: Dims = Dims([0, 0, 0, 0, 1, 0]);
const V: Dims = Dims([2, 1, -3, 0, -1, 0]);
const W: Dims = Dims([2, 1, -3, 0, 0, 0]);

fn key(node: &str, port: &str) -> PortKey {
    (node.into(), port.into())
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-8 * b.abs().max(1.0), "{a} != {b}");
}
fn controls() -> MaterialGraphSolve {
    MaterialGraphSolve {
        maximum_unknowns: 8,
        maximum_terms: 32,
        maximum_evaluations: 400,
        maximum_iterations: 24,
        equation_tolerance: 1e-10,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn with_cx<R>(run: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 81,
                kernel_id: 4,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &gate)
    })
}

fn hyperelastic_graph(mu_pa: f64) -> ConstitutiveGraph {
    let source = fs_blake3::hash_bytes(b"synthetic registered hyperelastic traction fixture");
    let card = ConstitutiveModelCard {
        law: LawId("neo-hookean-log".into()),
        law_version: 1,
        parameters: [
            (
                "mu".into(),
                LawParameter {
                    value: mu_pa,
                    dims: PA,
                },
            ),
            (
                "lambda".into(),
                LawParameter {
                    value: 0.0,
                    dims: PA,
                },
            ),
            (
                "stretch_limit".into(),
                LawParameter {
                    value: 5.0,
                    dims: Dims::NONE,
                },
            ),
        ]
        .into(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained(),
        sources: vec![source],
        provenance: Provenance {
            source: "synthetic traction fixture".into(),
            license: "CC0-1.0".into(),
            artifact: Some(source),
        },
    };
    let mut registry = LawRegistry::new();
    register_hyperelastic_laws(&mut registry);
    let mut graph = ConstitutiveGraph::new();
    graph
        .add_node("solid", registry.instantiate("solid", &card).unwrap())
        .unwrap();
    graph
}

#[test]
fn native_model_card_solves_prescribed_traction_and_material_only_swap() {
    with_cx(|cx, _| {
        let mut displacements = Vec::new();
        for mu in [3e5, 6e5] {
            let graph = hyperelastic_graph(mu);
            let fixed = (1..9)
                .map(|i| (key("solid", &format!("H_{}{}", i / 3, i % 3)), 0.0))
                .collect();
            let load = 1.5e5;
            let residual = MaterialGraphResidual::new(
                cx,
                &graph,
                graph.state_schema().version(),
                &[],
                fixed,
                vec![MaterialUnknown {
                    input: key("solid", "H_00"),
                    dims: Dims::NONE,
                    scale_si: 0.1,
                }],
                vec![MaterialEquation {
                    terms: vec![(1.0, MaterialValue::Output(key("solid", "P_00")))],
                    dims: PA,
                    target_si: load,
                    scale_si: load,
                }],
                controls(),
            )
            .unwrap();
            let solution = residual.solve(vec![0.0]).unwrap();
            // lambda=0 gives P=mu*(stretch-1/stretch), a quadratic independent
            // of the native constitutive implementation and Newton driver.
            let ratio = load / mu;
            let stretch = 0.5 * (ratio + (ratio * ratio + 4.0).sqrt());
            close(solution.external[&key("solid", "H_00")], stretch - 1.0);
            close(solution.graph.outputs[&key("solid", "P_00")], load);
            assert!(solution.graph.next_state.is_empty());
            assert!(
                solution.evaluations > 1 && solution.evaluations <= controls().maximum_evaluations
            );
            displacements.push(solution.external[&key("solid", "H_00")]);
        }
        assert!(displacements[1] < displacements[0]);
    });
}

fn declaration(
    name: &str,
    inputs: &[(&str, Dims)],
    outputs: &[(&str, Dims)],
    state: &[&str],
) -> NodeDeclaration {
    let ports = |values: &[(&str, Dims)]| {
        values
            .iter()
            .map(|(name, dims)| Port {
                name: (*name).into(),
                dims: *dims,
                parity: TimeParity::Even,
            })
            .collect()
    };
    NodeDeclaration {
        law: LawId(name.into()),
        law_version: 1,
        role: NodeRole::InternalMemory,
        inputs: ports(inputs),
        outputs: ports(outputs),
        state_slots: state.iter().map(|s| (*s).into()).collect(),
        state_schema_version: 1,
        calibration: ValidityDomain::unconstrained(),
        differentiability: Differentiability::Smooth,
        energy: EnergyBehavior::Empirical,
        tangent_claimed: true,
    }
}

struct Resistor {
    declaration: NodeDeclaration,
    reference_ohm: f64,
    alpha_per_k: f64,
}
impl Resistor {
    fn resistance(&self, t: f64) -> f64 {
        self.reference_ohm * (1.0 + self.alpha_per_k * (t - 300.0))
    }
}
impl LawNode for Resistor {
    fn declaration(&self) -> &NodeDeclaration {
        &self.declaration
    }
    fn evaluate(&self, _state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        if !(0.0..=500.0).contains(&inputs[0]) {
            return Err(GraphError::CardMismatch {
                law: "synthetic-resistor".into(),
                obligation: "temperature outside fixture source",
            });
        }
        let r = self.resistance(inputs[0]);
        let voltage = inputs[1] * r;
        let power = voltage * inputs[1];
        Ok(NodeOutput {
            outputs: vec![voltage, power],
            next_state: vec![],
            dissipation_rate: Some(power),
        })
    }
    fn tangent(&self, _state: &[f64], inputs: &[f64]) -> Option<Vec<f64>> {
        let dr = self.reference_ohm * self.alpha_per_k;
        let i = inputs[1];
        let r = self.resistance(inputs[0]);
        Some(vec![i * dr, r, i * i * dr, 2.0 * i * r])
    }
}
struct Caloric {
    declaration: NodeDeclaration,
    dt_over_capacity: f64,
}
impl LawNode for Caloric {
    fn declaration(&self) -> &NodeDeclaration {
        &self.declaration
    }
    fn evaluate(&self, state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        let next = state[0] + self.dt_over_capacity * inputs[0];
        Ok(NodeOutput {
            outputs: vec![next],
            next_state: vec![next],
            dissipation_rate: None,
        })
    }
    fn tangent(&self, _state: &[f64], _inputs: &[f64]) -> Option<Vec<f64>> {
        Some(vec![self.dt_over_capacity])
    }
}

fn electrothermal_graph() -> ConstitutiveGraph {
    let mut graph = ConstitutiveGraph::new();
    let mut resistor = declaration(
        "synthetic-resistor",
        &[("temperature", K), ("current", A)],
        &[("voltage", V), ("heat", W)],
        &[],
    );
    resistor.energy = EnergyBehavior::NonNegativeDissipation;
    graph
        .add_node(
            "resistor",
            Box::new(Resistor {
                declaration: resistor,
                reference_ohm: 2.0,
                alpha_per_k: 0.001,
            }),
        )
        .unwrap();
    graph
        .add_node(
            "thermal",
            Box::new(Caloric {
                declaration: declaration(
                    "finite-capacity",
                    &[("heat", W)],
                    &[("temperature", K)],
                    &["temperature"],
                ),
                dt_over_capacity: 0.2,
            }),
        )
        .unwrap();
    graph
        .connect(Edge {
            from: "resistor".into(),
            from_port: "heat".into(),
            to: "thermal".into(),
            to_port: "heat".into(),
        })
        .unwrap();
    graph
}
fn unknowns() -> Vec<MaterialUnknown> {
    vec![
        MaterialUnknown {
            input: key("resistor", "temperature"),
            dims: K,
            scale_si: 300.0,
        },
        MaterialUnknown {
            input: key("resistor", "current"),
            dims: A,
            scale_si: 5.0,
        },
    ]
}
fn equations() -> Vec<MaterialEquation> {
    vec![
        MaterialEquation {
            terms: vec![(1.0, MaterialValue::Output(key("resistor", "voltage")))],
            dims: V,
            target_si: 10.0,
            scale_si: 10.0,
        },
        MaterialEquation {
            terms: vec![
                (1.0, MaterialValue::Input(key("resistor", "temperature"))),
                (-1.0, MaterialValue::Output(key("thermal", "temperature"))),
            ],
            dims: K,
            target_si: 0.0,
            scale_si: 300.0,
        },
    ]
}

#[test]
fn closed_electrothermal_feedback_uses_one_committed_history_for_every_probe() {
    with_cx(|cx, _| {
        let graph = electrothermal_graph();
        let original = vec![300.0];
        let residual = MaterialGraphResidual::new(
            cx,
            &graph,
            graph.state_schema().version(),
            &original,
            BTreeMap::new(),
            unknowns(),
            equations(),
            controls(),
        )
        .unwrap();
        let mut first = [0.0; 2];
        let mut repeated = [0.0; 2];
        residual.residual(&[1.0, 1.0], &mut first);
        residual.residual(&[1.0, 1.0], &mut repeated);
        assert_eq!(first, repeated);
        let mut action = [0.0; 2];
        residual.jacobian_apply(&[1.0, 1.0], &[0.25, -0.5], &mut action);
        // Exact derivative of (I R(T)/10-1, (T-Told-.2 I²R(T))/300).
        close(action[0], 0.3 * 0.25 - 0.5);
        close(action[1], 0.99 * 0.25 - (1.0 / 15.0) * -0.5);
        assert_eq!(original, vec![300.0]);
        let solution = residual.solve(vec![1.0, 1.0]).unwrap();
        let rise = 20.0 / (1.0 + 1.04_f64.sqrt());
        close(
            solution.external[&key("resistor", "temperature")],
            300.0 + rise,
        );
        close(
            solution.external[&key("resistor", "current")],
            10.0 / (2.0 * (1.0 + 0.001 * rise)),
        );
        close(solution.graph.next_state[0], 300.0 + rise);
        close(solution.graph.total_dissipation * 0.2, rise);
        assert_eq!(original, vec![300.0]);
        // Continuing uses precisely the accepted memory, once, without any
        // probe count appearing in physical time or accumulated heat.
        let continued = MaterialGraphResidual::new(
            cx,
            &graph,
            graph.state_schema().version(),
            &solution.graph.next_state,
            BTreeMap::new(),
            unknowns(),
            equations(),
            controls(),
        )
        .unwrap()
        .solve(solution.coordinates.clone())
        .unwrap();
        close(
            continued.graph.next_state[0] - solution.graph.next_state[0],
            0.2 * continued.graph.total_dissipation,
        );
    });
}

#[test]
fn equation_dimensions_evaluation_budget_and_source_failure_cannot_become_convergence() {
    with_cx(|cx, gate| {
        let graph = electrothermal_graph();
        let state = [300.0];
        let version = graph.state_schema().version();
        let mut wrong = equations();
        wrong[0].dims = W;
        assert!(matches!(
            MaterialGraphResidual::new(
                cx,
                &graph,
                version,
                &state,
                BTreeMap::new(),
                unknowns(),
                wrong,
                controls()
            ),
            Err(MaterialGraphError::Invalid(_))
        ));
        let mut limited = controls();
        limited.maximum_evaluations = 1;
        let residual = MaterialGraphResidual::new(
            cx,
            &graph,
            version,
            &state,
            BTreeMap::new(),
            unknowns(),
            equations(),
            limited,
        )
        .unwrap();
        let mut output = [0.0; 2];
        residual.residual(&[1.0, 1.0], &mut output);
        assert!(output.iter().all(|x| x.is_finite()));
        residual.jacobian_apply(&[1.0, 1.0], &[1.0, 0.0], &mut output);
        assert!(output.iter().all(|x| x.is_nan()));
        assert!(matches!(
            residual.take_failure(),
            Some(MaterialGraphError::EvaluationBudget)
        ));
        assert!(residual.has_failed());
        residual.residual(&[1.0, 1.0], &mut output);
        assert!(output.iter().all(|x| x.is_nan()));
        assert!(residual.has_failed());
        let source = MaterialGraphResidual::new(
            cx,
            &graph,
            version,
            &state,
            BTreeMap::new(),
            unknowns(),
            equations(),
            controls(),
        )
        .unwrap();
        assert!(matches!(
            source.solve(vec![2.0, 1.0]),
            Err(MaterialGraphError::Graph(GraphError::CardMismatch { .. }))
        ));
        assert_eq!(state, [300.0]);
        gate.request();
        assert!(matches!(
            MaterialGraphResidual::new(
                cx,
                &graph,
                version,
                &state,
                BTreeMap::new(),
                unknowns(),
                equations(),
                controls()
            ),
            Err(MaterialGraphError::Cancelled)
        ));
    });
}
