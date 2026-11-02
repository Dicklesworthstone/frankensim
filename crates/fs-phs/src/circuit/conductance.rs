//! Grounded reciprocal multi-terminal conductance in the native circuit DAE.
use super::{CircuitError, CircuitGraph};
use crate::linear_exact::positive_heat::factor_resistance;
use std::sync::Arc;

/// Dense preparation and every query have at most 32 independent voltages.
pub const MAX_CONDUCTANCE_BLOCK_COORDINATES: usize = 32;
/// Maximum number of prepared blocks in one bounded assembly.
pub const MAX_CONDUCTANCE_BLOCKS: usize = 64;

/// Complete passive law in a caller-selected terminal order.
///
/// `nodes[0]` is the reference; `I[1..] = G (V[1..] - V[0])` and
/// `I[0] = -sum(I[1..])`. Positive current enters the conducting body.
/// The original reciprocal G [S] is retained, including either off-diagonal
/// sign. Its checked positive factor supplies independent dissipated power;
/// it does not replace G by a fitted resistor graph or change circuit storage.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductanceBlock {
    nodes: Arc<[usize]>,
    matrix_s: Arc<[f64]>,
    reference_s: Arc<[f64]>,
    root: Arc<[f64]>,
}

impl ConductanceBlock {
    /// Prepare an exactly reciprocal, finite, passive grounded matrix.
    /// Numerically unresolved negative modes refuse without clipping.
    pub fn try_new(nodes: Vec<usize>, matrix_s: Vec<f64>) -> Result<Self, CircuitError> {
        Self::try_new_polled(nodes, matrix_s, || Ok(()))
    }

    /// Same preparation with cooperative cancellation before and after the
    /// existing bounded resistance-factor kernel. That kernel is synchronous:
    /// at most 32 coordinates, cubic factor storage and quartic scalar work.
    /// Every later query uses quadratic work and linear temporary storage.
    pub fn try_new_polled(
        nodes: Vec<usize>,
        matrix_s: Vec<f64>,
        mut poll: impl FnMut() -> Result<(), CircuitError>,
    ) -> Result<Self, CircuitError> {
        poll()?;
        if !(2..=MAX_CONDUCTANCE_BLOCK_COORDINATES + 1).contains(&nodes.len())
            || nodes
                .iter()
                .enumerate()
                .any(|(i, node)| nodes[..i].contains(node))
        {
            return Err(invalid("conductance needs 2..=33 distinct terminal nodes"));
        }
        let n = nodes.len() - 1;
        if matrix_s.len() != n * n || matrix_s.iter().any(|value| !value.is_finite()) {
            return Err(invalid("complete finite grounded conductance required"));
        }
        for i in 0..n {
            for j in 0..i {
                if matrix_s[i * n + j] != matrix_s[j * n + i] {
                    return Err(invalid("conductance must be exactly reciprocal"));
                }
            }
        }
        let root = factor_resistance(&matrix_s, n)?;
        poll()?;
        // These are precisely the reference row and column of B G B^T.
        // Compensated sums preserve a reference which has exact zero current.
        let mut reference_s = vec![0.0; n + 1];
        for i in 0..n {
            reference_s[i + 1] = -sum(matrix_s[i * n..(i + 1) * n].iter().copied())?;
        }
        reference_s[0] = -sum(reference_s[1..].iter().copied())?;
        if reference_s[0] < 0.0
            || (reference_s[0] == 0.0 && reference_s[1..].iter().any(|v| *v != 0.0))
        {
            return Err(invalid("conductance reference loss is unresolved"));
        }
        Ok(Self {
            nodes: nodes.into(),
            matrix_s: matrix_s.into(),
            reference_s: reference_s.into(),
            root: root.into(),
        })
    }

    /// Original graph nodes; the first is this field's voltage reference.
    pub fn nodes(&self) -> &[usize] {
        &self.nodes
    }

    /// Original grounded G [S], row-major in non-reference terminal order.
    pub fn matrix_s(&self) -> &[f64] {
        &self.matrix_s
    }

    /// Complete terminal currents [A], with exact reference current closure
    /// to the represented non-reference current sum.
    pub fn terminal_currents_a(&self, voltage_v: &[f64]) -> Result<Vec<f64>, CircuitError> {
        let relative = self.relative_voltages(voltage_v)?;
        let n = relative.len();
        let mut currents = vec![0.0; n + 1];
        for (i, row) in self.matrix_s.chunks_exact(n).enumerate() {
            currents[i + 1] = sum(row.iter().zip(&relative).map(|(g, v)| g * v))?;
        }
        currents[0] = -sum(currents[1..].iter().copied())?;
        Ok(currents)
    }

    /// Independent physical loss [W] from the admitted positive squares.
    /// No heat is inferred from a circuit storage change or ledger residual.
    pub fn power_w(&self, voltage_v: &[f64]) -> Result<f64, CircuitError> {
        let relative = self.relative_voltages(voltage_v)?;
        let mut powers = Vec::with_capacity(relative.len());
        for row in self.root.chunks_exact(relative.len()) {
            let effort = sum(row.iter().zip(&relative).map(|(g, v)| g * v))?;
            powers.push(effort * effort);
        }
        sum(powers)
    }

    fn relative_voltages(&self, voltage_v: &[f64]) -> Result<Vec<f64>, CircuitError> {
        if voltage_v.len() != self.nodes.len() || voltage_v.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite terminal voltages required"));
        }
        voltage_v[1..]
            .iter()
            .map(|v| finite(v - voltage_v[0]))
            .collect()
    }

    fn coefficient(&self, i: usize, j: usize) -> f64 {
        match (i, j) {
            (0, _) => self.reference_s[j],
            (_, 0) => self.reference_s[i],
            _ => self.matrix_s[(i - 1) * (self.nodes.len() - 1) + j - 1],
        }
    }

    pub(super) fn admit_nodes(&self, incident: &mut [usize]) -> Result<(), CircuitError> {
        for (i, &node) in self.nodes.iter().enumerate() {
            if node >= incident.len() {
                return Err(invalid("conductance terminal node out of range"));
            }
            // A zero block or a genuinely absent reference must not conceal
            // the existing isolated-node refusal.
            if self.coefficient(i, i) > 0.0 {
                incident[node] = 1;
            }
        }
        Ok(())
    }

    pub(super) fn stamp(&self, r: &mut [f64], n: usize, off_v: usize) -> Result<(), CircuitError> {
        for (i, &a) in self.nodes.iter().enumerate() {
            if a == 0 {
                continue;
            }
            let row = off_v + a - 1;
            for (j, &b) in self.nodes[..=i].iter().enumerate() {
                if b == 0 {
                    continue;
                }
                let col = off_v + b - 1;
                let value = self.coefficient(i, j);
                r[row * n + col] = finite(r[row * n + col] + value)?;
                if row != col {
                    r[col * n + row] = finite(r[col * n + row] + value)?;
                }
            }
        }
        Ok(())
    }
}

/// Bound the new dense seam before entering the original native assembler.
/// The original assembler and integrator remain the sole Kirchhoff owners.
pub(super) fn admit_circuit(
    graph: &CircuitGraph,
    blocks: &[ConductanceBlock],
) -> Result<(), CircuitError> {
    if graph.node_count > 64
        || graph.branches.len() > 128
        || graph.transformers.len() > 32
        || blocks.len() > MAX_CONDUCTANCE_BLOCKS
    {
        return Err(invalid(
            "conductance circuit topology exceeds dense work bounds",
        ));
    }
    let mut states = graph.node_count.saturating_sub(1) + graph.transformers.len();
    for &(_, _, branch) in &graph.branches {
        match branch {
            super::Branch::Capacitor { .. } => states += 2,
            super::Branch::Inductor { .. } => states += 1,
            super::Branch::VoltageSource { port } | super::Branch::CurrentSource { port } => {
                if port >= 32 {
                    return Err(invalid("conductance circuit exceeds 32 source ports"));
                }
                states += usize::from(matches!(branch, super::Branch::VoltageSource { .. }));
            }
            super::Branch::Resistor { .. } => {}
        }
    }
    if states > 256 {
        return Err(invalid("conductance circuit exceeds 256 descriptor states"));
    }
    Ok(())
}

fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64, CircuitError> {
    let (mut total, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        finite(value)?;
        let next = finite(total + value)?;
        correction += if total.abs() >= value.abs() {
            (total - next) + value
        } else {
            (value - next) + total
        };
        total = next;
    }
    finite(total + correction)
}

fn finite(value: f64) -> Result<f64, CircuitError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("conductance arithmetic is unrepresentable"))
    }
}

fn invalid(what: &'static str) -> CircuitError {
    CircuitError::Invalid { what }
}
