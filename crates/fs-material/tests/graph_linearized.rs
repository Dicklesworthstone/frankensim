//! Chain rules must differentiate the same trial and committed internal history.
use fs_evidence::ValidityDomain;
use fs_matdb::LawId;
use fs_material::graph::*;
use fs_qty::Dims;
use std::collections::BTreeMap;

struct MemoryLaw {
    declaration: NodeDeclaration,
    bad: bool,
}
impl MemoryLaw {
    fn new(two_inputs: bool) -> Self {
        let port = |name: &str| Port {
            name: name.into(),
            dims: Dims::NONE,
            parity: TimeParity::Even,
        };
        Self {
            declaration: NodeDeclaration {
                law: LawId("synthetic-history-polynomial".into()),
                law_version: 1,
                role: NodeRole::InternalMemory,
                inputs: if two_inputs {
                    vec![port("x"), port("y")]
                } else {
                    vec![port("x")]
                },
                outputs: vec![port("value")],
                state_slots: vec!["history".into()],
                state_schema_version: 1,
                calibration: ValidityDomain::unconstrained(),
                differentiability: Differentiability::Smooth,
                energy: EnergyBehavior::Empirical,
                tangent_claimed: true,
            },
            bad: false,
        }
    }
}
impl LawNode for MemoryLaw {
    fn declaration(&self) -> &NodeDeclaration {
        &self.declaration
    }
    fn evaluate(&self, state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        Ok(NodeOutput {
            outputs: if self.bad {
                vec![]
            } else {
                vec![state[0] * inputs[0] * inputs[0] + 3.0 * inputs.get(1).copied().unwrap_or(0.0)]
            },
            next_state: vec![state[0] + inputs[0]],
            dissipation_rate: None,
        })
    }
    fn tangent(&self, state: &[f64], inputs: &[f64]) -> Option<Vec<f64>> {
        Some(if inputs.len() == 1 {
            vec![2.0 * state[0] * inputs[0]]
        } else {
            vec![2.0 * state[0] * inputs[0], 3.0]
        })
    }
}
fn key(node: &str, port: &str) -> (String, String) {
    (node.into(), port.into())
}
fn fixture() -> (
    ConstitutiveGraph,
    u64,
    Vec<f64>,
    BTreeMap<(String, String), f64>,
) {
    let mut graph = ConstitutiveGraph::new();
    // Reverse insertion order exercises a genuine topological chain.
    graph.add_node("b", Box::new(MemoryLaw::new(true))).unwrap();
    graph
        .add_node("a", Box::new(MemoryLaw::new(false)))
        .unwrap();
    graph
        .connect(Edge {
            from: "a".into(),
            from_port: "value".into(),
            to: "b".into(),
            to_port: "x".into(),
        })
        .unwrap();
    let (version, state) = graph.state_schema().encode(&[&[5.0], &[7.0]]).unwrap();
    (
        graph,
        version,
        state,
        BTreeMap::from([(key("a", "x"), 2.0), (key("b", "y"), 4.0)]),
    )
}

#[test]
fn graph_chain_rule_matches_analytic_and_independent_finite_difference() {
    let (graph, version, state, external) = fixture();
    let direction = BTreeMap::from([(key("a", "x"), 0.5), (key("b", "y"), -1.0)]);
    let linear = graph
        .linearize(version, &state, &external, &direction, || false)
        .unwrap();
    assert_eq!(
        linear.value,
        graph.execute(version, &state, &external).unwrap()
    );
    assert_eq!(linear.value.outputs[&key("b", "value")], 3932.0);
    assert_eq!(linear.output_direction[&key("b", "value")], 3917.0);
    assert_eq!(linear.value.next_state, vec![33.0, 9.0]);
    assert_eq!(state, vec![5.0, 7.0]);
    let epsilon = 1e-5;
    let sample = |sign: f64| {
        let shifted = external
            .iter()
            .map(|(key, value)| {
                (
                    key.clone(),
                    value + sign * epsilon * direction.get(key).copied().unwrap_or(0.0),
                )
            })
            .collect();
        graph.execute(version, &state, &shifted).unwrap().outputs[&key("b", "value")]
    };
    assert!(((sample(1.0) - sample(-1.0)) / (2.0 * epsilon) - 3917.0).abs() < 1e-6);
    assert_eq!(graph.external_input_dims("b", "x"), None);
    assert_eq!(graph.external_input_dims("b", "y"), Some(Dims::NONE));
}

#[test]
fn sparse_external_directions_hold_other_inputs_fixed() {
    let (graph, version, state, external) = fixture();
    let direction = BTreeMap::from([(key("b", "y"), 2.0)]);
    let linear = graph
        .linearize(version, &state, &external, &direction, || false)
        .unwrap();
    assert_eq!(linear.output_direction[&key("a", "value")], 0.0);
    assert_eq!(linear.output_direction[&key("b", "value")], 6.0);
    assert!(matches!(
        graph.linearize(
            version,
            &state,
            &external,
            &BTreeMap::from([(key("b", "x"), 1.0)]),
            || false
        ),
        Err(GraphError::UnknownEndpoint { .. })
    ));
}

#[test]
fn interruption_and_malformed_laws_return_no_partial_history() {
    let (graph, version, state, external) = fixture();
    let mut polls = 0;
    assert!(matches!(
        graph.linearize(version, &state, &external, &BTreeMap::new(), || {
            polls += 1;
            polls >= 5
        }),
        Err(GraphError::Cancelled)
    ));
    assert_eq!(state, vec![5.0, 7.0]);
    let mut broken = ConstitutiveGraph::new();
    let mut law = MemoryLaw::new(false);
    law.bad = true;
    broken.add_node("bad", Box::new(law)).unwrap();
    let (version, state) = broken.state_schema().encode(&[&[1.0]]).unwrap();
    assert!(matches!(
        broken.execute(version, &state, &BTreeMap::from([(key("bad", "x"), 1.0)])),
        Err(GraphError::ArityMismatch {
            which: "outputs",
            ..
        })
    ));
}

#[test]
fn nonfinite_state_input_and_tangent_actions_refuse() {
    let (graph, version, state, mut external) = fixture();
    external.insert(key("a", "x"), f64::NAN);
    assert!(matches!(
        graph.execute(version, &state, &external),
        Err(GraphError::NonFinite { .. })
    ));
    external.insert(key("a", "x"), 2.0);
    assert!(matches!(
        graph.linearize(
            version,
            &state,
            &external,
            &BTreeMap::from([(key("a", "x"), f64::MAX)]),
            || false
        ),
        Err(GraphError::NonFinite { .. })
    ));
}
