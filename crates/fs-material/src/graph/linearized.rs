//! The same DAG evaluation supplies values and matrix-free chain-rule actions.
use super::*;

/// One proposed constitutive update and its directional output derivative.
/// Internal history is held at the supplied committed state throughout. This
/// is the algorithmic input tangent, not differentiation through earlier steps.
#[derive(Debug, Clone, PartialEq)]
pub struct GraphLinearization {
    /// Original law values, proposed history and reported dissipation.
    pub value: GraphOutput,
    /// `d output / d external * direction`, keyed like `value.outputs`.
    pub output_direction: BTreeMap<(String, String), f64>,
}

impl ConstitutiveGraph {
    /// Dimensions of an input which has no internal driver. This is the legal
    /// boundary for an external solver unknown; internal edges stay DAG-only.
    #[must_use]
    pub fn external_input_dims(&self, node: &str, port: &str) -> Option<Dims> {
        if self
            .edges
            .iter()
            .any(|edge| edge.to == node && edge.to_port == port)
        {
            return None;
        }
        self.nodes
            .get(node)?
            .declaration()
            .inputs
            .iter()
            .find(|p| p.name == port)
            .map(|p| p.dims)
    }

    /// Dimensions of a declared output, for typed external balance equations.
    #[must_use]
    pub fn output_dims(&self, node: &str, port: &str) -> Option<Dims> {
        self.nodes
            .get(node)?
            .declaration()
            .outputs
            .iter()
            .find(|p| p.name == port)
            .map(|p| p.dims)
    }

    /// Evaluate with checkpoints before/after every bounded law callback.
    /// A refusal returns no partial history. A law must bound its own work;
    /// these checkpoints cannot preempt an individual callback.
    pub fn execute_cancellable(
        &self,
        state_version: u64,
        state: &[f64],
        external: &BTreeMap<(String, String), f64>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GraphOutput, GraphError> {
        Ok(self
            .evaluate_direction(state_version, state, external, None, cancelled)?
            .value)
    }

    /// Propagate a direction through declared consistent tangents at exactly
    /// the same trial inputs and committed history as the value evaluation.
    /// No global Jacobian or finite-difference stencil is constructed. Omitted
    /// external direction entries are zero. Undeclared/driven directions,
    /// missing tangents, and malformed/nonfinite node results refuse.
    pub fn linearize(
        &self,
        state_version: u64,
        state: &[f64],
        external: &BTreeMap<(String, String), f64>,
        direction: &BTreeMap<(String, String), f64>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GraphLinearization, GraphError> {
        for ((node, port), &value) in direction {
            if self.external_input_dims(node, port).is_none() {
                return Err(GraphError::UnknownEndpoint {
                    reference: format!("external {node}.{port}"),
                });
            }
            finite(node, "input direction", &[value])?;
        }
        self.evaluate_direction(state_version, state, external, Some(direction), cancelled)
    }

    fn evaluate_direction(
        &self,
        state_version: u64,
        state: &[f64],
        external: &BTreeMap<(String, String), f64>,
        direction: Option<&BTreeMap<(String, String), f64>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GraphLinearization, GraphError> {
        let mut check = || {
            if cancelled() {
                Err(GraphError::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        finite("graph", "committed state", state)?;
        let schema = self.state_schema();
        let mut node_states = schema.decode(state_version, state)?;
        let order = self.topological_order()?;
        check()?;
        let mut produced = BTreeMap::new();
        let mut output_direction = BTreeMap::new();
        let mut total_dissipation = 0.0;
        for node_id in &order {
            check()?;
            let node = &self.nodes[node_id];
            let declaration = node.declaration();
            let mut inputs = Vec::with_capacity(declaration.inputs.len());
            let mut input_direction = Vec::new();
            for port in &declaration.inputs {
                let edge = self
                    .edges
                    .iter()
                    .find(|e| &e.to == node_id && e.to_port == port.name);
                let key = (node_id.clone(), port.name.clone());
                let value = edge
                    .map(|e| produced[&(e.from.clone(), e.from_port.clone())])
                    .or_else(|| external.get(&key).copied())
                    .ok_or_else(|| GraphError::UnfedPort {
                        node: node_id.clone(),
                        port: port.name.clone(),
                    })?;
                inputs.push(value);
                if let Some(d) = direction {
                    input_direction.push(edge.map_or_else(
                        || d.get(&key).copied().unwrap_or(0.0),
                        |e| output_direction[&(e.from.clone(), e.from_port.clone())],
                    ));
                }
            }
            finite(node_id, "inputs", &inputs)?;
            let slot = self
                .order
                .iter()
                .position(|id| id == node_id)
                .expect("admitted node");
            let result = node.evaluate(&node_states[slot], &inputs)?;
            check()?;
            arity(
                node_id,
                "outputs",
                declaration.outputs.len(),
                result.outputs.len(),
            )?;
            arity(
                node_id,
                "next state",
                declaration.state_slots.len(),
                result.next_state.len(),
            )?;
            finite(node_id, "outputs", &result.outputs)?;
            finite(node_id, "next state", &result.next_state)?;
            if direction.is_some() {
                let tangent_error = || GraphError::TangentClaimUnmet {
                    node: node_id.clone(),
                    law: declaration.law.0.clone(),
                };
                if !declaration.tangent_claimed {
                    return Err(tangent_error());
                }
                let tangent = node
                    .tangent(&node_states[slot], &inputs)
                    .ok_or_else(tangent_error)?;
                check()?;
                let count = declaration
                    .outputs
                    .len()
                    .checked_mul(inputs.len())
                    .ok_or_else(tangent_error)?;
                if tangent.len() != count {
                    return Err(tangent_error());
                }
                finite(node_id, "tangent", &tangent)?;
                for (port, row) in declaration
                    .outputs
                    .iter()
                    .zip(tangent.chunks_exact(inputs.len()))
                {
                    let value: f64 = row.iter().zip(&input_direction).map(|(a, b)| a * b).sum();
                    finite(node_id, "output direction", &[value])?;
                    output_direction.insert((node_id.clone(), port.name.clone()), value);
                }
            }
            if declaration.energy == EnergyBehavior::ReversibleTransfer
                && result.dissipation_rate != Some(0.0)
            {
                return Err(GraphError::CardMismatch {
                    law: declaration.law.0.clone(),
                    obligation: "reversible transfer must report exactly zero dissipation",
                });
            }
            if let Some(rate) = result.dissipation_rate {
                finite(node_id, "dissipation rate", &[rate])?;
                if matches!(
                    declaration.energy,
                    EnergyBehavior::NonNegativeDissipation | EnergyBehavior::StorageAndDissipation
                ) && rate < 0.0
                {
                    return Err(GraphError::NegativeDissipation {
                        node: node_id.clone(),
                        rate,
                    });
                }
                total_dissipation += rate;
                finite(node_id, "total dissipation", &[total_dissipation])?;
            }
            for (port, value) in declaration.outputs.iter().zip(result.outputs) {
                produced.insert((node_id.clone(), port.name.clone()), value);
            }
            node_states[slot] = result.next_state;
        }
        check()?;
        let per_node: Vec<&[f64]> = node_states.iter().map(Vec::as_slice).collect();
        let (_, next_state) = schema.encode(&per_node)?;
        Ok(GraphLinearization {
            value: GraphOutput {
                outputs: produced,
                next_state,
                total_dissipation,
            },
            output_direction,
        })
    }
}

fn finite(node: &str, field: &'static str, values: &[f64]) -> Result<(), GraphError> {
    if let Some(value) = values.iter().find(|v| !v.is_finite()) {
        return Err(GraphError::NonFinite {
            node: node.into(),
            field,
            bits: value.to_bits(),
        });
    }
    Ok(())
}
fn arity(node: &str, which: &'static str, expected: usize, found: usize) -> Result<(), GraphError> {
    if expected != found {
        return Err(GraphError::ArityMismatch {
            node: node.into(),
            which,
            expected,
            found,
        });
    }
    Ok(())
}
