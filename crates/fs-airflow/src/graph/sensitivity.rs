//! Implicit tangents and adjoints of a converged passive loss graph.
//!
//! For `dp = A q + B q |q|`, `g = 1/(A + 2 B |q|)` and a common
//! coefficient scale gives `dq = g d(dp) - g dp d(log scale)`. This reduces
//! to the original `q/2` response for quadratic edges and to `q` for Darcy
//! edges. Eliminating fixed pressures gives the symmetric
//! positive-definite, grounded weighted graph Laplacian. We solve that system,
//! not a differentiation of the primal iteration history. Storage is O(V+E).
//!
//! These are local derivatives of the declared nominal passive model, not
//! interval certificates, uncertainty propagation, or coupled thermal adjoints.
//! Fan-driven derivatives are available through `solve_with_fan_linearized`.
//! Zero-flow quadratic branches are explicitly refused. A positive linear
//! coefficient supplies its true finite zero-flow slope without regularization.
//! Mixed-edge coefficient directions scale both coefficients together; separate
//! linear/quadratic parameter directions are not represented by this API.
//! With gravity, `dp` here is the dissipative drive (static dp plus gravity
//! head). These derivatives hold gravity, positions and density samples fixed;
//! scaling resistance never scales the conservative gravitational head.

use core::fmt;
use std::cell::RefCell;

use fs_exec::Cx;
use fs_math::det;
use fs_qty::{Pressure, Qty, VolumetricFlowRate};
use fs_solver::{CheckedCgConfig, CheckedCgError, LinearOp, SolveReport, checked_cg};
use fs_sparse::precond::Precond;

use super::{
    FixedPressure, GraphError, GraphSolution, GraphSolveConfig, LossGraph, Sum,
    checkpoint, finite, periodic_checkpoint,
};

mod fan;
pub use fan::{GraphFanGradient, GraphFanLinearization, GraphFanSensitivityError, GraphFanTangent};

/// Derivative of a dimensionless objective with respect to pressure (1/Pa).
pub type PressureWeight = Qty<1, -1, 2, 0, 0, 0>;
/// Derivative of a dimensionless objective with respect to volume flow (s/m^3).
pub type FlowWeight = Qty<-3, 0, 1, 0, 0, 0>;

/// Explicit admission and work limits for the local derivative model.
#[derive(Debug, Clone, Copy)]
pub struct GraphSensitivityConfig {
    /// Refuse branches at or below this nonnegative absolute flow (m^3/s).
    /// Zero also admits stagnant branches with a positive linear coefficient.
    pub minimum_absolute_flow: VolumetricFlowRate,
    /// Shared iteration budget and recomputed-residual target for each solve.
    pub linear: CheckedCgConfig,
}

impl GraphSensitivityConfig {
    fn validate(self) -> Result<(), GraphSensitivityError> {
        let minimum = self.minimum_absolute_flow.value();
        if !(minimum.is_finite() && minimum >= 0.0) {
            return Err(GraphSensitivityError::InvalidInput("minimum flow must be finite and nonnegative"));
        }
        if self.linear.max_iterations == 0
            || !(self.linear.tolerance.is_finite()
                && self.linear.tolerance > 0.0 && self.linear.tolerance < 1.0)
        {
            return Err(GraphSensitivityError::InvalidInput("positive iteration budget and tolerance in (0, 1) required"));
        }
        Ok(())
    }
}

/// Input, primal, differentiability, interruption, or linear-solve refusal.
#[derive(Debug)]
pub enum GraphSensitivityError {
    /// A dimension, value, or derivative budget is invalid.
    InvalidInput(&'static str),
    /// Primal validation, convergence, arithmetic, or cancellation failed.
    Graph(GraphError),
    /// The pressure-to-flow law is not admitted for differentiation here.
    NonDifferentiableBranch {
        /// Original branch index.
        branch: usize,
        /// Signed nominal flow at the refused operating point.
        flow: VolumetricFlowRate,
    },
    /// The recomputed linear residual did not meet the requested tolerance.
    LinearDidNotConverge(SolveReport),
}

impl From<GraphError> for GraphSensitivityError {
    fn from(error: GraphError) -> Self { Self::Graph(error) }
}

impl fmt::Display for GraphSensitivityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "airflow graph sensitivity: {self:?}")
    }
}

impl std::error::Error for GraphSensitivityError {}

/// Directional derivative with respect to a dimensionless caller parameter.
#[derive(Debug, Clone)]
pub struct GraphTangent {
    /// Pressure derivatives in node order, including prescribed pressures.
    pub pressures: Vec<Pressure>,
    /// Signed flow derivatives in original branch order.
    pub flows: Vec<VolumetricFlowRate>,
    /// Net directional outflows in node order; fixed nodes exchange with reservoirs.
    pub node_outflows: Vec<VolumetricFlowRate>,
    /// Observed maximum free-node differentiated-conservation residual (m^3/s).
    pub max_node_imbalance: VolumetricFlowRate,
    /// Recomputed, floating-point linear residual; not a forward-error bound.
    pub linear_report: SolveReport,
}

/// Gradient of a dimensionless scalar objective, from one transposed solve.
#[derive(Debug, Clone)]
pub struct GraphGradient {
    /// Derivative for a common multiplicative resistance scale on each edge.
    /// For a pure law this is `dJ/d(log R)`; on mixed edges both A and B scale.
    pub log_resistances: Vec<f64>,
    /// `dJ/dp_boundary` in the original boundary-input order (1/Pa).
    pub boundary_pressures: Vec<PressureWeight>,
    /// Recomputed, floating-point adjoint residual; not a gradient-error bound.
    pub linear_report: SolveReport,
}

/// Immutable solved state and admitted linearization, reusable across objectives.
///
/// Constructed only by [`LossGraph::solve_linearized`]; mutable or externally
/// assembled `GraphSolution` values cannot bypass primal acceptance.
#[derive(Debug)]
pub struct GraphLinearization<'g> {
    graph: &'g LossGraph,
    solution: GraphSolution,
    free_nodes: Vec<usize>,
    free_index: Vec<Option<usize>>,
    conductance: Vec<f64>,
    resistance_response: Vec<f64>,
    inverse_diagonal: Vec<f64>,
    config: GraphSensitivityConfig,
}

impl LossGraph {
    /// Solve and admit implicit derivatives of the nominal passive graph.
    ///
    /// Each component is grounded by the primal boundary validation. Positive
    /// finite conductances then make the free-node Laplacian SPD. Admission
    /// refuses non-differentiable quadratic stagnation, caller-excluded low
    /// flows, and unrepresentable coefficients.
    pub fn solve_linearized(
        &self,
        boundaries: &[FixedPressure],
        primal: GraphSolveConfig,
        sensitivity: GraphSensitivityConfig,
        cx: &Cx<'_>,
    ) -> Result<GraphLinearization<'_>, GraphSensitivityError> {
        checkpoint(cx)?;
        sensitivity.validate()?;
        let solution = self.solve(boundaries, primal, cx)?;
        self.linearize_accepted(solution, sensitivity, cx)
    }

    // Private admission seam: callers here own the accepted primal solve. Never
    // expose this for externally assembled or mutated GraphSolution values.
    fn linearize_accepted(
        &self, solution: GraphSolution, sensitivity: GraphSensitivityConfig, cx: &Cx<'_>,
    ) -> Result<GraphLinearization<'_>, GraphSensitivityError> {
        checkpoint(cx)?;
        sensitivity.validate()?;
        let boundaries = &solution.boundaries;
        let mut fixed = vec![false; self.node_count()];
        for (index, boundary) in boundaries.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            fixed[boundary.node] = true;
        }
        let mut free_nodes = Vec::new();
        let mut free_index = vec![None; self.node_count()];
        for (node, &is_fixed) in fixed.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            if !is_fixed {
                free_index[node] = Some(free_nodes.len());
                free_nodes.push(node);
            }
        }
        let mut conductance = Vec::with_capacity(self.branches.len());
        let mut resistance_response = Vec::with_capacity(self.branches.len());
        let mut diagonal = vec![Sum::default(); free_nodes.len()];
        for (edge, branch) in solution.branches.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            let q = branch.flow.value();
            let linear = self.linear_resistance[edge];
            if q.abs() <= sensitivity.minimum_absolute_flow.value()
                && (linear == 0.0 || sensitivity.minimum_absolute_flow.value() > 0.0)
            {
                return Err(GraphSensitivityError::NonDifferentiableBranch {
                    branch: edge, flow: branch.flow,
                });
            }
            let g = if linear == 0.0 {
                // Preserve the original separate-root quadratic evaluation.
                finite((0.5 / self.sqrt_resistance[edge])
                    / det::sqrt(branch.dissipative_pressure_drop.value().abs()), "branch conductance")?
            } else {
                let quadratic = finite(branch.loss.resistance.value() * q.abs(),
                    "quadratic derivative coefficient")?;
                let scale = linear.max(quadratic);
                finite((1.0 / (linear / scale + 2.0 * (quadratic / scale))) / scale,
                    "mixed branch conductance")?
            };
            let response = if linear == 0.0 { 0.5 * q }
                else { finite(g * branch.dissipative_pressure_drop.value(), "resistance scale response")? };
            if g <= 0.0 || (q != 0.0 && response == 0.0) {
                return Err(GraphError::NonFinite("branch derivative underflowed").into());
            }
            conductance.push(g);
            resistance_response.push(response);
            for node in [branch.from, branch.to] {
                if let Some(index) = free_index[node] { diagonal[index].add(g)?; }
            }
        }
        let mut inverse_diagonal = Vec::with_capacity(free_nodes.len());
        for (index, sum) in diagonal.into_iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            let inverse = finite(1.0 / sum.value()?, "inverse Laplacian diagonal")?;
            if inverse <= 0.0 {
                return Err(GraphError::NonFinite("inverse diagonal underflowed").into());
            }
            inverse_diagonal.push(inverse);
        }
        checkpoint(cx)?;
        Ok(GraphLinearization {
            graph: self, solution, free_nodes, free_index, conductance,
            resistance_response, inverse_diagonal, config: sensitivity,
        })
    }
}

impl GraphLinearization<'_> {
    /// Accepted primal result with its original physical metadata unchanged.
    #[must_use]
    pub fn solution(&self) -> &GraphSolution { &self.solution }

    /// Jacobian-vector product for logarithmic resistances and fixed pressures.
    ///
    /// Each direction scales all nonzero resistance coefficients of that edge:
    /// `dA/dt = A direction`, `dB/dt = B direction`. On a Darcy-only edge,
    /// negative log-permeability direction gives the log-resistance direction
    /// when viscosity and geometry are fixed. Boundary directions retain their
    /// original input order; the parameter `t` is dimensionless.
    pub fn tangent(
        &self,
        log_resistance_direction: &[f64],
        boundary_pressure_direction: &[Pressure],
        cx: &Cx<'_>,
    ) -> Result<GraphTangent, GraphSensitivityError> {
        checkpoint(cx)?;
        if log_resistance_direction.len() != self.graph.branches.len()
            || boundary_pressure_direction.len() != self.solution.boundaries.len()
        {
            return Err(GraphSensitivityError::InvalidInput("tangent direction length mismatch"));
        }
        let mut pressures = vec![0.0; self.graph.node_count()];
        for (index, (boundary, direction)) in self.solution.boundaries.iter()
            .zip(boundary_pressure_direction).enumerate()
        {
            periodic_checkpoint(index, cx)?;
            pressures[boundary.node] = finite(direction.value(), "boundary pressure direction")?;
        }
        let mut rhs = vec![Sum::default(); self.free_nodes.len()];
        for (edge, branch) in self.solution.branches.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            let direction = finite(log_resistance_direction[edge], "log-resistance direction")?;
            let imposed = finite(self.conductance[edge]
                * (pressures[branch.from] - pressures[branch.to])
                - self.resistance_response[edge] * direction, "imposed tangent flow")?;
            if let Some(index) = self.free_index[branch.from] { rhs[index].add(-imposed)?; }
            if let Some(index) = self.free_index[branch.to] { rhs[index].add(imposed)?; }
        }
        let (free_pressures, linear_report) = self.solve_rhs(rhs, cx)?;
        for (index, &node) in self.free_nodes.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            pressures[node] = free_pressures[index];
        }
        let mut flows = Vec::with_capacity(self.graph.branches.len());
        let mut outflows = vec![Sum::default(); self.graph.node_count()];
        for (edge, branch) in self.solution.branches.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            let flow = finite(self.conductance[edge]
                * (pressures[branch.from] - pressures[branch.to])
                - self.resistance_response[edge] * log_resistance_direction[edge], "tangent flow")?;
            flows.push(VolumetricFlowRate::new(flow));
            outflows[branch.from].add(flow)?;
            outflows[branch.to].add(-flow)?;
        }
        let mut node_outflows = Vec::with_capacity(outflows.len());
        let mut max_node_imbalance: f64 = 0.0;
        for (node, sum) in outflows.into_iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            let value = sum.value()?;
            if self.free_index[node].is_some() {
                max_node_imbalance = max_node_imbalance.max(value.abs());
            }
            node_outflows.push(VolumetricFlowRate::new(value));
        }
        let mut typed_pressures = Vec::with_capacity(pressures.len());
        for (node, value) in pressures.into_iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            typed_pressures.push(Pressure::new(value));
        }
        checkpoint(cx)?;
        Ok(GraphTangent {
            pressures: typed_pressures, flows, node_outflows,
            max_node_imbalance: VolumetricFlowRate::new(max_node_imbalance), linear_report,
        })
    }

    /// Reverse-mode gradient for `dJ = sum(w_p dp) + sum(w_q dq)`.
    ///
    /// Pressure weights are in node order, flow weights in branch order. This
    /// includes weights on prescribed pressures. Explicit objective dependence
    /// on R or other inputs must be added by the caller after this chain rule.
    pub fn adjoint(
        &self,
        pressure_weights: &[PressureWeight],
        flow_weights: &[FlowWeight],
        cx: &Cx<'_>,
    ) -> Result<GraphGradient, GraphSensitivityError> {
        checkpoint(cx)?;
        if pressure_weights.len() != self.graph.node_count()
            || flow_weights.len() != self.graph.branches.len()
        {
            return Err(GraphSensitivityError::InvalidInput("objective weight length mismatch"));
        }
        let mut rhs = vec![Sum::default(); self.graph.node_count()];
        for (node, weight) in pressure_weights.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            rhs[node].add(finite(weight.value(), "pressure objective weight")?)?;
        }
        for (edge, branch) in self.graph.branches.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            let weight = finite(flow_weights[edge].value(), "flow objective weight")?;
            let value = finite(self.conductance[edge] * weight, "adjoint right hand side")?;
            rhs[branch.from].add(value)?;
            rhs[branch.to].add(-value)?;
        }
        let mut free_rhs = Vec::with_capacity(self.free_nodes.len());
        for (index, &node) in self.free_nodes.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            free_rhs.push(rhs[node]);
        }
        let (lambda_free, linear_report) = self.solve_rhs(free_rhs, cx)?;
        let mut lambda = vec![0.0; self.graph.node_count()];
        for (index, &node) in self.free_nodes.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            lambda[node] = lambda_free[index];
        }
        let mut log_resistances = Vec::with_capacity(self.graph.branches.len());
        let mut boundary_gradient = vec![Sum::default(); self.graph.node_count()];
        for (node, weight) in pressure_weights.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            boundary_gradient[node].add(weight.value())?;
        }
        for (edge, branch) in self.solution.branches.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            let corrected = finite(flow_weights[edge].value()
                - (lambda[branch.from] - lambda[branch.to]), "adjoint edge weight")?;
            log_resistances.push(finite(-self.resistance_response[edge] * corrected,
                "log-resistance gradient")?);
            let value = finite(self.conductance[edge] * corrected, "boundary pressure gradient")?;
            boundary_gradient[branch.from].add(value)?;
            boundary_gradient[branch.to].add(-value)?;
        }
        let mut boundary_pressures = Vec::with_capacity(self.solution.boundaries.len());
        for (index, boundary) in self.solution.boundaries.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            boundary_pressures.push(PressureWeight::new(boundary_gradient[boundary.node].value()?));
        }
        checkpoint(cx)?;
        Ok(GraphGradient { log_resistances, boundary_pressures, linear_report })
    }

    fn solve_rhs(
        &self, sums: Vec<Sum>, cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, SolveReport), GraphSensitivityError> {
        let mut rhs = Vec::with_capacity(sums.len());
        for (index, sum) in sums.into_iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            rhs.push(sum.value()?);
        }
        let operator = LinearizedOperator { state: self, cx, failure: RefCell::new(None) };
        let solved = checked_cg(&operator, &operator, &rhs, self.config.linear,
            |_| operator.check()).map_err(|error| match error {
                CheckedCgError::InvalidInput(message) => GraphSensitivityError::InvalidInput(message),
                CheckedCgError::Interrupted(error) => GraphSensitivityError::Graph(error),
            })?;
        if !solved.report.converged_euclidean() {
            return Err(GraphSensitivityError::LinearDidNotConverge(solved.report));
        }
        checkpoint(cx)?;
        Ok((solved.x, solved.report))
    }
}

// The solver's infallible operator/preconditioner traits need a refusal latch.
// Checked CG polls it after operator work; NaNs prevent accepting partial work
// inside an inner batch. All operator traversal itself checkpoints every tile.
struct LinearizedOperator<'a, 'g, 'cx> {
    state: &'a GraphLinearization<'g>,
    cx: &'a Cx<'cx>,
    failure: RefCell<Option<GraphError>>,
}

impl LinearizedOperator<'_, '_, '_> {
    fn check(&self) -> Result<(), GraphError> {
        if let Some(error) = self.failure.borrow_mut().take() { return Err(error); }
        checkpoint(self.cx)
    }

    fn publish(&self, result: Result<(), GraphError>, output: &mut [f64]) {
        if let Err(error) = result {
            let mut failure = self.failure.borrow_mut();
            if failure.is_none() { *failure = Some(error); }
            output.fill(f64::NAN);
        }
    }

    fn apply_laplacian(&self, x: &[f64], y: &mut [f64]) -> Result<(), GraphError> {
        checkpoint(self.cx)?;
        for (index, &node) in self.state.free_nodes.iter().enumerate() {
            periodic_checkpoint(index, self.cx)?;
            let mut sum = Sum::default();
            for (ordinal, &edge) in self.state.graph.adjacency[node].iter().enumerate() {
                periodic_checkpoint(ordinal, self.cx)?;
                let other = self.state.graph.other(edge, node);
                let neighbor = self.state.free_index[other].map_or(0.0, |j| x[j]);
                sum.add(self.state.conductance[edge] * (x[index] - neighbor))?;
            }
            y[index] = sum.value()?;
        }
        Ok(())
    }

    fn apply_diagonal(&self, r: &[f64], z: &mut [f64]) -> Result<(), GraphError> {
        checkpoint(self.cx)?;
        for (index, &inverse) in self.state.inverse_diagonal.iter().enumerate() {
            periodic_checkpoint(index, self.cx)?;
            z[index] = finite(r[index] * inverse, "Jacobi preconditioner")?;
        }
        Ok(())
    }
}

impl LinearOp for LinearizedOperator<'_, '_, '_> {
    fn n(&self) -> usize { self.state.free_nodes.len() }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.publish(self.apply_laplacian(x, y), y);
    }
}

impl Precond for LinearizedOperator<'_, '_, '_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        self.publish(self.apply_diagonal(r, z), z);
    }
}

#[cfg(test)]
mod tests;
