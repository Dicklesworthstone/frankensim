//! Steady passive quadratic-loss networks with arbitrary graph topology.
//!
//! An edge directed from `u` to `v` obeys `p[u] - p[v] = R Q |Q|`.
//! The direction is bookkeeping: negative flow is permitted. Free nodes obey
//! volume conservation; every connected component needs a prescribed pressure.
//! Node-wise bracketed equilibration avoids dividing by a zero-flow derivative.
//!
//! Results are nominal numerical estimates, NOT interval certificates. Retained
//! coefficient uncertainty and validity cards are not propagated or validated by
//! this solver. It supplies no CFD, thermal-feedback, transient, adjoint, or
//! cross-ISA accuracy claim. Unlike `EnclosureNetwork`, this generic graph does
//! not infer or require a leakage branch: callers must supply every physical
//! path. Work uses fixed node/edge order, `fs_math::det`, explicit iteration
//! limits, and checkpoints in each node solve and bounded edge traversals.

use core::fmt;
use std::collections::BTreeSet;

use fs_exec::Cx;
use fs_math::det;
use fs_qty::{Pressure, VolumetricFlowRate};

use crate::LossElement;

mod fan;
pub use fan::{GraphFanConfig, GraphFanError, GraphOperatingPoint};

/// An oriented passive branch; its orientation does not constrain flow sign.
#[derive(Debug, Clone, PartialEq)]
pub struct GraphBranch {
    /// Upstream index for positive flow.
    pub from: usize,
    /// Downstream index for positive flow.
    pub to: usize,
    /// Retained nominal coefficient, uncertainty, source, and validity card.
    pub loss: LossElement,
}

/// A prescribed nodal pressure in Pa, relative to any common reference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FixedPressure {
    /// Node index.
    pub node: usize,
    /// Fixed pressure relative to the caller's reference.
    pub pressure: Pressure,
}

/// Explicit iteration and volume-conservation budgets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GraphSolveConfig {
    /// Maximum complete sweeps over all free nodes.
    pub max_sweeps: usize,
    /// Maximum bracket bisections in each node update.
    pub max_node_iterations: usize,
    /// Absolute conservation tolerance in m^3/s, strictly positive.
    pub absolute_flow_tolerance: VolumetricFlowRate,
    /// Relative conservation tolerance in [0, 1).
    pub relative_flow_tolerance: f64,
}

impl GraphSolveConfig {
    fn validate(self) -> Result<(), GraphError> {
        if self.max_sweeps == 0 || self.max_node_iterations == 0 {
            return Err(GraphError::InvalidInput("iteration budgets must be positive"));
        }
        let absolute = self.absolute_flow_tolerance.value();
        if !(absolute.is_finite() && absolute > 0.0) {
            return Err(GraphError::InvalidInput("absolute flow tolerance must be finite and positive"));
        }
        if !(self.relative_flow_tolerance.is_finite()
            && (0.0..1.0).contains(&self.relative_flow_tolerance))
        {
            return Err(GraphError::InvalidInput("relative flow tolerance must be in [0, 1)"));
        }
        Ok(())
    }

    fn threshold(self, magnitude: f64) -> Result<f64, GraphError> {
        finite(
            self.absolute_flow_tolerance.value()
                + self.relative_flow_tolerance * (0.5 * magnitude),
            "conservation tolerance",
        )
    }
}

/// A solved edge, retaining the coefficient's physical authority unchanged.
#[derive(Debug, Clone, PartialEq)]
pub struct GraphBranchFlow {
    /// Positive-flow origin.
    pub from: usize,
    /// Positive-flow destination.
    pub to: usize,
    /// Original coefficient and its unpromoted physical metadata.
    pub loss: LossElement,
    /// Signed `pressure[from] - pressure[to]` in Pa.
    pub pressure_drop: Pressure,
    /// Signed nominal flow in m^3/s.
    pub flow: VolumetricFlowRate,
}

/// Residual-bearing nominal solution, not a certified enclosure.
#[derive(Debug, Clone, PartialEq)]
pub struct GraphSolution {
    /// Pressures in node-index order.
    pub pressures: Vec<Pressure>,
    /// Branches in input order, including zero and reverse flows.
    pub branches: Vec<GraphBranchFlow>,
    /// Net flow leaving each node. Reservoir nodes supply or absorb this flow.
    pub node_outflows: Vec<VolumetricFlowRate>,
    /// Original prescribed pressures.
    pub boundaries: Vec<FixedPressure>,
    /// Completed node sweeps; zero means the initial state already balanced.
    pub sweeps: usize,
    /// Maximum absolute free-node conservation residual in m^3/s.
    pub max_node_imbalance: VolumetricFlowRate,
    /// Signed sum of the reservoir exchanges in m^3/s.
    pub boundary_imbalance: VolumetricFlowRate,
}

/// Structured input, arithmetic, interruption, and convergence failures.
#[derive(Debug)]
pub enum GraphError {
    /// Invalid input/configuration; the message names the violated condition.
    InvalidInput(&'static str),
    /// Invalid edge endpoint, self-loop, coefficient, or identity.
    InvalidBranch {
        /// Index in the input edge sequence.
        branch: usize,
        /// Violated condition.
        reason: &'static str,
    },
    /// A component has no fixed pressure and therefore no pressure reference.
    UnreferencedNode(usize),
    /// A context checkpoint refused further work.
    Interrupted,
    /// Floating-point arithmetic could not represent the requested operation.
    NonFinite(&'static str),
    /// The iteration budget ended without meeting conservation tolerances.
    DidNotConverge {
        /// Completed sweeps.
        sweeps: usize,
        /// Maximum free-node imbalance in m^3/s.
        max_node_imbalance: VolumetricFlowRate,
        /// Signed net boundary exchange in m^3/s.
        boundary_imbalance: VolumetricFlowRate,
    },
}

impl fmt::Display for GraphError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "airflow graph: {self:?}")
    }
}

impl std::error::Error for GraphError {}

/// Validated reciprocal quadratic-loss topology. Parallel edges are allowed.
#[derive(Debug, Clone)]
pub struct LossGraph {
    branches: Vec<GraphBranch>,
    adjacency: Vec<Vec<usize>>,
    sqrt_resistance: Vec<f64>,
}

impl LossGraph {
    /// Validate topology and retain the original loss data without promotion.
    pub fn new(node_count: usize, branches: Vec<GraphBranch>) -> Result<Self, GraphError> {
        if node_count == 0 {
            return Err(GraphError::InvalidInput("node count must be positive"));
        }
        let mut names = BTreeSet::new();
        let mut sqrt_resistance = Vec::with_capacity(branches.len());
        for (index, branch) in branches.iter().enumerate() {
            let invalid = |reason| GraphError::InvalidBranch { branch: index, reason };
            if branch.from >= node_count || branch.to >= node_count || branch.from == branch.to {
                return Err(invalid("endpoints must be distinct valid node indices"));
            }
            if branch.loss.name.trim().is_empty() || !names.insert(branch.loss.name.clone()) {
                return Err(invalid("branch names must be nonblank and unique"));
            }
            let resistance = branch.loss.resistance.value();
            if !(resistance.is_finite() && resistance > 0.0) {
                return Err(invalid("resistance must be finite and positive"));
            }
            let uncertainty = branch.loss.uncertainty_rel;
            if !(uncertainty.is_finite() && (0.0..1.0).contains(&uncertainty)) {
                return Err(invalid("resistance uncertainty must be in [0, 1)"));
            }
            let root = det::sqrt(resistance);
            if !(root.is_finite() && root > 0.0) {
                return Err(invalid("resistance square root is not representable"));
            }
            sqrt_resistance.push(root);
        }
        let mut adjacency = vec![Vec::new(); node_count];
        for (index, branch) in branches.iter().enumerate() {
            adjacency[branch.from].push(index);
            adjacency[branch.to].push(index);
        }
        Ok(Self { branches, adjacency, sqrt_resistance })
    }

    /// Number of nodes, including explicit isolated nodes.
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.adjacency.len()
    }

    /// Immutable input edges, preserving all physical metadata.
    #[must_use]
    pub fn branches(&self) -> &[GraphBranch] {
        &self.branches
    }

    /// Solve arbitrary looped or disconnected-but-anchored pressure networks.
    ///
    /// Each free node and the summed reservoir exchange must satisfy
    /// `|imbalance| <= absolute + relative * (sum of absolute flows)/2`.
    /// Exhaustion, disconnected pressure references, and interruption refuse;
    /// an unconverged field is never returned as a successful solution.
    pub fn solve(
        &self,
        boundaries: &[FixedPressure],
        config: GraphSolveConfig,
        cx: &Cx<'_>,
    ) -> Result<GraphSolution, GraphError> {
        checkpoint(cx)?;
        config.validate()?;
        let fixed = self.fixed_pressures(boundaries, cx)?;
        let low = boundaries.iter().map(|b| b.pressure.value()).fold(f64::INFINITY, f64::min);
        let high = boundaries.iter().map(|b| b.pressure.value()).fold(f64::NEG_INFINITY, f64::max);
        let midpoint = 0.5 * low + 0.5 * high;
        let mut pressures: Vec<f64> = fixed.iter().map(|value| value.unwrap_or(midpoint)).collect();
        let mut state = FlowState::new(self.node_count(), self.branches.len());
        for sweep in 0..=config.max_sweeps {
            if self.evaluate(&pressures, &fixed, config, cx, &mut state)? {
                return self.finish(pressures, boundaries, state, sweep, cx);
            }
            if sweep == config.max_sweeps {
                return Err(state.nonconvergence(sweep));
            }
            self.relax_flat_components(&mut pressures, &fixed, config.max_node_iterations, cx)?;
            for node in 0..self.node_count() {
                checkpoint(cx)?;
                if fixed[node].is_none() {
                    pressures[node] = self.equilibrate(node, &pressures, config.max_node_iterations, cx)?;
                }
            }
            self.relax_flat_components(&mut pressures, &fixed, config.max_node_iterations, cx)?;
        }
        unreachable!("inclusive sweep range always returns at its limit")
    }

    fn fixed_pressures(
        &self, boundaries: &[FixedPressure], cx: &Cx<'_>,
    ) -> Result<Vec<Option<f64>>, GraphError> {
        if boundaries.is_empty() {
            return Err(GraphError::InvalidInput("at least one prescribed pressure is required"));
        }
        let mut fixed = vec![None; self.node_count()];
        let mut visited = vec![false; self.node_count()];
        let mut queue = Vec::new();
        for boundary in boundaries {
            checkpoint(cx)?;
            if boundary.node >= self.node_count() || !boundary.pressure.value().is_finite() {
                return Err(GraphError::InvalidInput("boundary node or pressure is invalid"));
            }
            if fixed[boundary.node].replace(boundary.pressure.value()).is_some() {
                return Err(GraphError::InvalidInput("a node has duplicate pressure boundaries"));
            }
            visited[boundary.node] = true;
            queue.push(boundary.node);
        }
        let mut cursor = 0;
        while cursor < queue.len() {
            checkpoint(cx)?;
            let node = queue[cursor];
            cursor += 1;
            for (ordinal, &edge) in self.adjacency[node].iter().enumerate() {
                periodic_checkpoint(ordinal, cx)?;
                let other = self.other(edge, node);
                if !visited[other] {
                    visited[other] = true;
                    queue.push(other);
                }
            }
        }
        if let Some(node) = visited.iter().position(|seen| !seen) {
            return Err(GraphError::UnreferencedNode(node));
        }
        Ok(fixed)
    }

    fn other(&self, edge: usize, node: usize) -> usize {
        let branch = &self.branches[edge];
        if branch.from == node { branch.to } else { branch.from }
    }

    fn flow(&self, edge: usize, from_pressure: f64, to_pressure: f64) -> Result<f64, GraphError> {
        let difference = finite(from_pressure - to_pressure, "pressure difference")?;
        if difference == 0.0 {
            return Ok(0.0);
        }
        // Taking the roots separately also admits dp/R outside the f64 range
        // when its square root is still representable.
        let magnitude = finite(det::sqrt(difference.abs()) / self.sqrt_resistance[edge], "branch flow")?;
        if magnitude == 0.0 {
            return Err(GraphError::NonFinite("nonzero branch flow underflowed"));
        }
        Ok(magnitude.copysign(difference))
    }

    fn node_balance(
        &self, node: usize, pressure: f64, pressures: &[f64], cx: &Cx<'_>,
    ) -> Result<f64, GraphError> {
        let mut sum = Sum::default();
        for (ordinal, &edge) in self.adjacency[node].iter().enumerate() {
            periodic_checkpoint(ordinal, cx)?;
            sum.add(self.flow(edge, pressure, pressures[self.other(edge, node)])?)?;
        }
        sum.value()
    }

    fn equilibrate(
        &self, node: usize, pressures: &[f64], iterations: usize, cx: &Cx<'_>,
    ) -> Result<f64, GraphError> {
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for (ordinal, &edge) in self.adjacency[node].iter().enumerate() {
            periodic_checkpoint(ordinal, cx)?;
            let value = pressures[self.other(edge, node)];
            low = low.min(value);
            high = high.max(value);
        }
        if low == high {
            return Ok(low);
        }
        let mut best = pressures[node];
        let mut best_residual = self.node_balance(node, best, pressures, cx)?.abs();
        for candidate in [low, high] {
            let residual = self.node_balance(node, candidate, pressures, cx)?.abs();
            if residual < best_residual {
                best = candidate;
                best_residual = residual;
            }
        }
        for _ in 0..iterations {
            checkpoint(cx)?;
            let middle = 0.5 * low + 0.5 * high;
            if middle <= low || middle >= high {
                break;
            }
            let residual = self.node_balance(node, middle, pressures, cx)?;
            if residual.abs() < best_residual {
                best = middle;
                best_residual = residual.abs();
            }
            if residual == 0.0 {
                return Ok(middle);
            }
            if residual > 0.0 { high = middle; } else { low = middle; }
        }
        Ok(best)
    }

    // Common shifts of near-equal-pressure components remove the zero-flow
    // stiffness that can stall point relaxation (for example on a dead leg).
    // Internal pressure differences are not smoothed. The shift minimizes the
    // same convex loss potential: its derivative is the component's net flow.
    fn relax_flat_components(
        &self, pressures: &mut [f64], fixed: &[Option<f64>], iterations: usize,
        cx: &Cx<'_>,
    ) -> Result<(), GraphError> {
        let low_pressure = pressures.iter().copied().fold(f64::INFINITY, f64::min);
        let high_pressure = pressures.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let scale = det::sqrt(f64::EPSILON);
        let flat_limit = high_pressure * scale - low_pressure * scale;
        let mut parents: Vec<usize> = (0..self.node_count()).collect();
        for (edge, branch) in self.branches.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            if fixed[branch.from].is_none() && fixed[branch.to].is_none()
                && (pressures[branch.from] - pressures[branch.to]).abs() <= flat_limit
            {
                let a = representative(&mut parents, branch.from);
                let b = representative(&mut parents, branch.to);
                parents[a.max(b)] = a.min(b);
            }
        }
        let mut groups = vec![Vec::new(); self.node_count()];
        for (node, boundary) in fixed.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            let root = representative(&mut parents, node);
            if boundary.is_none() { groups[root].push(node); }
        }
        for nodes in groups.iter().filter(|nodes| nodes.len() > 1) {
            checkpoint(cx)?;
            let root = parents[nodes[0]];
            let mut cut = Vec::new();
            let mut low = f64::INFINITY;
            let mut high = f64::NEG_INFINITY;
            for &node in nodes {
                checkpoint(cx)?;
                for (ordinal, &edge) in self.adjacency[node].iter().enumerate() {
                    periodic_checkpoint(ordinal, cx)?;
                    let other = self.other(edge, node);
                    if parents[other] != root {
                        cut.push((edge, node, other));
                        let difference = finite(pressures[other] - pressures[node], "block pressure bracket")?;
                        low = low.min(difference);
                        high = high.max(difference);
                    }
                }
            }
            if cut.is_empty() { return Err(GraphError::UnreferencedNode(nodes[0])); }
            let mut best = 0.0;
            let mut best_residual = self.block_balance(&cut, best, pressures, cx)?.abs();
            for candidate in [low, high] {
                let residual = self.block_balance(&cut, candidate, pressures, cx)?.abs();
                if residual < best_residual { best = candidate; best_residual = residual; }
            }
            for _ in 0..iterations {
                checkpoint(cx)?;
                if best_residual == 0.0 { break; }
                let middle = 0.5 * low + 0.5 * high;
                if middle <= low || middle >= high { break; }
                let residual = self.block_balance(&cut, middle, pressures, cx)?;
                if residual.abs() < best_residual { best = middle; best_residual = residual.abs(); }
                if residual > 0.0 { high = middle; } else { low = middle; }
            }
            for (ordinal, &node) in nodes.iter().enumerate() {
                periodic_checkpoint(ordinal, cx)?;
                pressures[node] = finite(pressures[node] + best, "block pressure update")?;
            }
        }
        Ok(())
    }

    fn block_balance(
        &self, cut: &[(usize, usize, usize)], shift: f64, pressures: &[f64], cx: &Cx<'_>,
    ) -> Result<f64, GraphError> {
        let mut sum = Sum::default();
        for (ordinal, &(edge, inside, outside)) in cut.iter().enumerate() {
            periodic_checkpoint(ordinal, cx)?;
            let shifted = finite(pressures[inside] + shift, "block trial pressure")?;
            sum.add(self.flow(edge, shifted, pressures[outside])?)?;
        }
        sum.value()
    }

    fn evaluate(
        &self, pressures: &[f64], fixed: &[Option<f64>], config: GraphSolveConfig,
        cx: &Cx<'_>, state: &mut FlowState,
    ) -> Result<bool, GraphError> {
        state.sums.fill(Sum::default());
        state.magnitudes.fill(Sum::default());
        state.max_node_imbalance = 0.0;
        for (edge, branch) in self.branches.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            let flow = self.flow(edge, pressures[branch.from], pressures[branch.to])?;
            state.flows[edge] = flow;
            state.sums[branch.from].add(flow)?;
            state.sums[branch.to].add(-flow)?;
            state.magnitudes[branch.from].add(flow.abs())?;
            state.magnitudes[branch.to].add(flow.abs())?;
        }
        let mut converged = true;
        let mut boundary_sum = Sum::default();
        let mut boundary_magnitude = Sum::default();
        for (node, boundary) in fixed.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            let residual = state.sums[node].value()?;
            state.balances[node] = residual;
            if boundary.is_some() {
                boundary_sum.add(residual)?;
                boundary_magnitude.add(residual.abs())?;
            } else {
                state.max_node_imbalance = state.max_node_imbalance.max(residual.abs());
                converged &= residual.abs() <= config.threshold(state.magnitudes[node].value()?)?;
            }
        }
        state.boundary_imbalance = boundary_sum.value()?;
        converged &= state.boundary_imbalance.abs() <= config.threshold(boundary_magnitude.value()?)?;
        checkpoint(cx)?;
        Ok(converged)
    }

    fn finish(
        &self, pressures: Vec<f64>, boundaries: &[FixedPressure], state: FlowState,
        sweeps: usize, cx: &Cx<'_>,
    ) -> Result<GraphSolution, GraphError> {
        let mut branches = Vec::with_capacity(self.branches.len());
        for (edge, branch) in self.branches.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            branches.push(GraphBranchFlow {
                from: branch.from,
                to: branch.to,
                loss: branch.loss.clone(),
                pressure_drop: Pressure::new(pressures[branch.from] - pressures[branch.to]),
                flow: VolumetricFlowRate::new(state.flows[edge]),
            });
        }
        let solution = GraphSolution {
            pressures: pressures.into_iter().map(Pressure::new).collect(),
            branches,
            node_outflows: state.balances.into_iter().map(VolumetricFlowRate::new).collect(),
            boundaries: boundaries.to_vec(),
            sweeps,
            max_node_imbalance: VolumetricFlowRate::new(state.max_node_imbalance),
            boundary_imbalance: VolumetricFlowRate::new(state.boundary_imbalance),
        };
        checkpoint(cx)?;
        Ok(solution)
    }
}

fn representative(parents: &mut [usize], node: usize) -> usize {
    let mut root = node;
    while parents[root] != root { root = parents[root]; }
    let mut current = node;
    while parents[current] != current {
        let next = parents[current];
        parents[current] = root;
        current = next;
    }
    root
}

fn finite(value: f64, stage: &'static str) -> Result<f64, GraphError> {
    if value.is_finite() { Ok(value) } else { Err(GraphError::NonFinite(stage)) }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), GraphError> {
    cx.checkpoint().map_err(|_| GraphError::Interrupted)
}

fn periodic_checkpoint(index: usize, cx: &Cx<'_>) -> Result<(), GraphError> {
    if index % 64 == 0 { checkpoint(cx)?; }
    Ok(())
}

#[derive(Debug, Clone, Copy, Default)]
struct Sum { sum: f64, correction: f64 }

impl Sum {
    fn add(&mut self, value: f64) -> Result<(), GraphError> {
        let next = finite(self.sum + value, "flow accumulation")?;
        let correction = if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.correction = finite(self.correction + correction, "compensated flow accumulation")?;
        self.sum = next;
        Ok(())
    }

    fn value(self) -> Result<f64, GraphError> {
        finite(self.sum + self.correction, "compensated flow sum")
    }
}

struct FlowState {
    flows: Vec<f64>,
    sums: Vec<Sum>,
    magnitudes: Vec<Sum>,
    balances: Vec<f64>,
    max_node_imbalance: f64,
    boundary_imbalance: f64,
}

impl FlowState {
    fn new(nodes: usize, branches: usize) -> Self {
        Self {
            flows: vec![0.0; branches],
            sums: vec![Sum::default(); nodes],
            magnitudes: vec![Sum::default(); nodes],
            balances: vec![0.0; nodes],
            max_node_imbalance: 0.0,
            boundary_imbalance: 0.0,
        }
    }

    fn nonconvergence(&self, sweeps: usize) -> GraphError {
        GraphError::DidNotConverge {
            sweeps,
            max_node_imbalance: VolumetricFlowRate::new(self.max_node_imbalance),
            boundary_imbalance: VolumetricFlowRate::new(self.boundary_imbalance),
        }
    }
}

#[cfg(test)]
mod tests;
