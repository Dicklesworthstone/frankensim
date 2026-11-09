//! Implicit derivatives of the original, converged Carreau pipe network.
//!
//! For each accepted radial profile let Q be signed volume flow, p its
//! constitutively reconstructed pressure drop and G=dQ/dp. Radial momentum
//! gives Q=pi R^3 H(pR/(2L)), hence dQ/dln(R)=3Q+pG and dQ/dln(L)=-pG.
//! With s=dln(a_T)/dT, the ORIGINAL Vertical shift gives dQ/dT=-s pG;
//! VerticalAndRate instead gives -s Q. The latter follows because
//! tau=a_T mu_ref(a_T gamma) gamma=tau_ref(a_T gamma).
//! The circular-pipe WRMS identity is given in Sochi (2015), equations 9-11:
//! <https://arxiv.org/pdf/1502.03314>. Shift conventions are those retained by
//! the material owner, not interchangeable fitted temperature laws.
//!
//! These identities reuse the accepted profile and the original shift slope;
//! there is no frozen apparent viscosity, radial replay, finite difference,
//! material refit or differentiation of Newton/continuation iterations.
//! Each tangent or adjoint uses the parent's grounded conductance Laplacian
//! and checked CG. The retained coefficients and storage are O(nodes+pipes).
//!
//! This is a local derivative of the declared nominal steady model. Source
//! support must contain a neighborhood of the accepted state, including the
//! retained physical rate ceiling under temperature shifts. A small linear
//! residual does not bound primal, quadrature, model or gradient error. No
//! transient/thermal-evolution adjoint or finite design-step admission is implied.

use fs_material::fluid_transport::rheology::ViscosityShiftMode;

use super::{
    BoundaryLayout, CarreauGraphControls, CarreauGraphError, CarreauGraphSolution,
    CarreauPipeGraph, CheckedCgConfig, CheckedCgError, ConductanceLaplacian, Cx, FixedPressure,
    FlowInjection, GraphError, Pressure, RefCell, SolveReport, Sum, VolumetricFlowRate, checked_cg,
    checkpoint, finite, periodic_checkpoint,
};
use crate::graph::sensitivity::{FlowWeight, PressureWeight};

/// Input derivatives with respect to one dimensionless caller parameter.
/// Geometry and temperatures follow original branch order; prescribed ports
/// follow their original input order, including repeated nodal injections.
#[derive(Clone, Copy, Debug)]
pub struct CarreauGraphDirection<'a> {
    /// d ln(R)/dt, with all other pipe declarations held fixed.
    pub log_radii: &'a [f64],
    /// d ln(L)/dt, with all other pipe declarations held fixed.
    pub log_lengths: &'a [f64],
    /// dT/dt [K], selecting the ORIGINAL temperature-shift convention.
    pub temperatures_k: &'a [f64],
    pub boundary_pressures: &'a [Pressure],
    pub injections: &'a [VolumetricFlowRate],
}

/// Accepted directional state; no failed linear iterate is returned.
#[derive(Clone, Debug)]
pub struct CarreauGraphTangent {
    pub pressures: Vec<Pressure>,
    pub flows: Vec<VolumetricFlowRate>,
    pub node_outflows: Vec<VolumetricFlowRate>,
    /// Net reservoir-flow derivative in original pressure-boundary order.
    /// A prescribed injection at a fixed node subtracts here without changing
    /// the pressure-controlled pipe solution.
    pub reservoir_outflows: Vec<VolumetricFlowRate>,
    /// Maximum free-node |d(outflow)-d(injection)| [m³/s].
    pub max_node_imbalance: VolumetricFlowRate,
    pub linear_report: SolveReport,
}

/// Gradient of a dimensionless objective, using one transposed solve.
/// The objective supplies node-pressure and branch-flow weights. Any direct
/// dependence on geometry, temperature or authored ports must be added by the
/// caller; this result is the chain rule through the hydraulic state only.
#[derive(Clone, Debug)]
pub struct CarreauGraphGradient {
    pub log_radii: Vec<f64>,
    pub log_lengths: Vec<f64>,
    /// dJ/dT for each original branch, in 1/K.
    pub temperatures_per_k: Vec<f64>,
    /// dJ/dp in original pressure-boundary order, in 1/Pa.
    pub boundary_pressures: Vec<PressureWeight>,
    /// dJ/dQ in original injection-input order, in s/m³. Repeated ports at one
    /// free node share a derivative; fixed-node ports have zero hydraulic
    /// derivative because their pressure reservoir absorbs the extra input.
    pub injections: Vec<FlowWeight>,
    pub linear_report: SolveReport,
}

/// Immutable accepted solution and original-law derivatives, reusable across
/// directions/objectives. Only a successful original primal solve constructs
/// this owner; mutable/external CarreauGraphSolution values cannot admit it.
#[derive(Debug)]
pub struct CarreauGraphLinearization<'g> {
    graph: &'g CarreauPipeGraph,
    solution: CarreauGraphSolution,
    layout: BoundaryLayout,
    reference_boundary: Vec<usize>,
    conductance: Vec<f64>,
    // [dQ/dln(R), dQ/dln(L), dQ/dT], all from the SAME accepted profile.
    responses: Vec<[f64; 3]>,
    inverse_diagonal: Vec<f64>,
    linear: CheckedCgConfig,
}

impl CarreauPipeGraph {
    pub fn solve_linearized(
        &self,
        boundaries: &[FixedPressure],
        primal: CarreauGraphControls,
        linear: CheckedCgConfig,
        cx: &Cx<'_>,
    ) -> Result<CarreauGraphLinearization<'_>, CarreauGraphError> {
        self.solve_with_injections_linearized(boundaries, &[], primal, linear, cx)
    }

    /// Solve and retain local sensitivities of the original radial source law.
    /// Rest and reversed flow use their true finite differential conductance.
    /// A source/physical rate boundary or temperature-domain endpoint refuses
    /// a two-sided derivative instead of extending material support.
    pub fn solve_with_injections_linearized(
        &self,
        boundaries: &[FixedPressure],
        injections: &[FlowInjection],
        primal: CarreauGraphControls,
        linear: CheckedCgConfig,
        cx: &Cx<'_>,
    ) -> Result<CarreauGraphLinearization<'_>, CarreauGraphError> {
        checkpoint(cx)?;
        if linear.max_iterations == 0
            || !linear.tolerance.is_finite()
            || linear.tolerance <= 0.0
            || linear.tolerance >= 1.0
        {
            return Err(GraphError::InvalidInput(
                "Carreau sensitivity requires finite linear tolerance in (0,1) and work",
            )
            .into());
        }
        let solution = self.solve_with_injections(boundaries, injections, primal, cx)?;
        let layout = self.boundary_layout(boundaries, injections, cx)?;
        let reference_boundary = component_references(self, boundaries, cx)?;
        let mut conductance = Vec::with_capacity(self.branches.len());
        let mut responses = Vec::with_capacity(self.branches.len());
        for (branch, accepted) in solution.branches.iter().enumerate() {
            periodic_checkpoint(branch, cx)?;
            let pipe = accepted.profile.pipe();
            let flow = accepted.profile.flow();
            let source = |source| CarreauGraphError::Pipe { branch, source };
            let invalid = |reason| GraphError::InvalidBranch { branch, reason };
            let temperature = pipe.temperature_k();
            let [low, high] = pipe.law().temperature_shift().validity_k();
            if temperature <= low
                || temperature >= high
                || flow.wall_shear_rate_per_s.abs() >= pipe.maximum_laminar_wall_rate_per_s()
            {
                return Err(invalid(
                    "Carreau sensitivity needs interior temperature and wall-rate support",
                )
                .into());
            }
            let original = pipe.law().evaluate(temperature, 0.0).map_err(source)?;
            // Preparing a perturbed pipe must retain its WHOLE authored rate
            // envelope, not only the current operating shear rate.
            if pipe.law().mode() == ViscosityShiftMode::VerticalAndRate {
                let ceiling =
                    pipe.maximum_laminar_wall_rate_per_s() * original.temperature().factor();
                if !ceiling.is_finite() || ceiling >= pipe.law().curve().rate_support_per_s()[1] {
                    return Err(invalid(
                        "Carreau temperature derivative needs interior reduced-rate ceiling",
                    )
                    .into());
                }
            }
            let g = finite(
                accepted.differential_conductance_m3_s_pa,
                "Carreau sensitivity conductance",
            )?;
            if g <= 0.0 {
                return Err(
                    invalid("Carreau differential conductance must remain positive").into(),
                );
            }
            let pressure_response = product(g, flow.pressure_drop_pa)?;
            let temperature_response = product(
                -original.temperature().log_derivative_per_k(),
                match pipe.law().mode() {
                    ViscosityShiftMode::Vertical => pressure_response,
                    ViscosityShiftMode::VerticalAndRate => flow.volumetric_flow_m3_s,
                },
            )?;
            conductance.push(g);
            responses.push([
                finite(
                    3.0 * flow.volumetric_flow_m3_s + pressure_response,
                    "Carreau radius response",
                )?,
                -pressure_response,
                temperature_response,
            ]);
        }
        let mut inverse_diagonal = Vec::with_capacity(layout.free_nodes.len());
        for (index, &node) in layout.free_nodes.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            let mut sum = Sum::default();
            for (ordinal, &edge) in self.adjacency[node].iter().enumerate() {
                periodic_checkpoint(ordinal, cx)?;
                sum.add(conductance[edge])?;
            }
            let inverse = finite(1.0 / sum.value()?, "Carreau sensitivity diagonal")?;
            if inverse <= 0.0 {
                return Err(GraphError::NonFinite("Carreau sensitivity diagonal underflow").into());
            }
            inverse_diagonal.push(inverse);
        }
        checkpoint(cx)?;
        Ok(CarreauGraphLinearization {
            graph: self,
            solution,
            layout,
            reference_boundary,
            conductance,
            responses,
            inverse_diagonal,
            linear,
        })
    }
}

// Cache the first authored pressure port in each already-grounded component.
// This traversal runs once when admitting the immutable linearization.
fn component_references(
    graph: &CarreauPipeGraph,
    boundaries: &[FixedPressure],
    cx: &Cx<'_>,
) -> Result<Vec<usize>, GraphError> {
    let mut reference = vec![usize::MAX; graph.node_count()];
    let mut queue = Vec::new();
    for (index, boundary) in boundaries.iter().enumerate() {
        periodic_checkpoint(index, cx)?;
        if reference[boundary.node] != usize::MAX {
            continue;
        }
        queue.clear();
        queue.push(boundary.node);
        reference[boundary.node] = index;
        let mut cursor = 0;
        while cursor < queue.len() {
            periodic_checkpoint(cursor, cx)?;
            let node = queue[cursor];
            for (ordinal, &edge) in graph.adjacency[node].iter().enumerate() {
                periodic_checkpoint(ordinal, cx)?;
                let other = graph.other(edge, node);
                if reference[other] == usize::MAX {
                    reference[other] = index;
                    queue.push(other);
                }
            }
            cursor += 1;
        }
    }
    Ok(reference)
}

fn product(a: f64, b: f64) -> Result<f64, GraphError> {
    let value = finite(a * b, "Carreau constitutive parameter derivative")?;
    if value == 0.0 && a != 0.0 && b != 0.0 {
        return Err(GraphError::NonFinite(
            "Carreau parameter derivative underflow",
        ));
    }
    Ok(value)
}

impl CarreauGraphLinearization<'_> {
    pub const fn solution(&self) -> &CarreauGraphSolution {
        &self.solution
    }

    /// Solve B_free*dQ = d(injection), with the original radial tangent and
    /// exact pipe-parameter responses. O(nodes+pipes) operator work/storage.
    pub fn tangent(
        &self,
        direction: CarreauGraphDirection<'_>,
        cx: &Cx<'_>,
    ) -> Result<CarreauGraphTangent, CarreauGraphError> {
        checkpoint(cx)?;
        let edges = self.graph.branches.len();
        if direction.log_radii.len() != edges
            || direction.log_lengths.len() != edges
            || direction.temperatures_k.len() != edges
            || direction.boundary_pressures.len() != self.solution.boundaries.len()
            || direction.injections.len() != self.solution.injections.len()
        {
            return Err(
                GraphError::InvalidInput("Carreau tangent direction length mismatch").into(),
            );
        }
        let mut pressures = vec![0.0; self.graph.node_count()];
        for (index, (boundary, value)) in self
            .solution
            .boundaries
            .iter()
            .zip(direction.boundary_pressures)
            .enumerate()
        {
            periodic_checkpoint(index, cx)?;
            let value = finite(value.value(), "Carreau boundary direction")?;
            let reference = finite(
                direction.boundary_pressures[self.reference_boundary[boundary.node]].value(),
                "Carreau component pressure direction",
            )?;
            pressures[boundary.node] = finite(value - reference, "centered pressure direction")?;
        }
        let mut injection = vec![Sum::default(); self.graph.node_count()];
        for (index, (port, value)) in self
            .solution
            .injections
            .iter()
            .zip(direction.injections)
            .enumerate()
        {
            periodic_checkpoint(index, cx)?;
            injection[port.node].add(finite(value.value(), "Carreau injection direction")?)?;
        }
        let mut rhs = injection.clone();
        let mut intrinsic = Vec::with_capacity(edges);
        for (edge, branch) in self.graph.branches.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            let mut response = Sum::default();
            for (coefficient, value) in self.responses[edge].into_iter().zip([
                direction.log_radii[edge],
                direction.log_lengths[edge],
                direction.temperatures_k[edge],
            ]) {
                response.add(product(
                    coefficient,
                    finite(value, "Carreau material/geometry direction")?,
                )?)?;
            }
            let response = response.value()?;
            intrinsic.push(response);
            let imposed = finite(
                self.conductance[edge] * (pressures[branch.from] - pressures[branch.to]) + response,
                "Carreau imposed tangent flow",
            )?;
            rhs[branch.from].add(-imposed)?;
            rhs[branch.to].add(imposed)?;
        }
        let (free_pressures, linear_report) = self.solve_rhs(&rhs, cx)?;
        for (index, &node) in self.layout.free_nodes.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            pressures[node] = free_pressures[index];
        }
        let mut flows = Vec::with_capacity(edges);
        let mut outflows = vec![Sum::default(); self.graph.node_count()];
        for (edge, branch) in self.graph.branches.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            let flow = finite(
                self.conductance[edge] * (pressures[branch.from] - pressures[branch.to])
                    + intrinsic[edge],
                "Carreau tangent flow",
            )?;
            flows.push(VolumetricFlowRate::new(flow));
            outflows[branch.from].add(flow)?;
            outflows[branch.to].add(-flow)?;
        }
        let mut node_outflows = Vec::with_capacity(outflows.len());
        let mut residuals = Vec::with_capacity(outflows.len());
        let mut max_node_imbalance: f64 = 0.0;
        for (node, sum) in outflows.into_iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            let outflow = sum.value()?;
            let residual = finite(
                outflow - injection[node].value()?,
                "Carreau differentiated continuity",
            )?;
            if !self.layout.fixed[node] {
                max_node_imbalance = max_node_imbalance.max(residual.abs());
            }
            node_outflows.push(VolumetricFlowRate::new(outflow));
            residuals.push(residual);
        }
        let mut reservoir_outflows = Vec::with_capacity(self.solution.boundaries.len());
        for (index, boundary) in self.solution.boundaries.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            reservoir_outflows.push(VolumetricFlowRate::new(residuals[boundary.node]));
        }
        // Flows use the centered solution, so a uniform component gauge
        // cannot consume the accuracy of a weak physical pressure direction.
        let mut physical_pressures = Vec::with_capacity(pressures.len());
        for (node, pressure) in pressures.into_iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            physical_pressures.push(Pressure::new(finite(
                pressure + direction.boundary_pressures[self.reference_boundary[node]].value(),
                "restored Carreau pressure direction",
            )?));
        }
        checkpoint(cx)?;
        Ok(CarreauGraphTangent {
            pressures: physical_pressures,
            flows,
            node_outflows,
            reservoir_outflows,
            max_node_imbalance: VolumetricFlowRate::new(max_node_imbalance),
            linear_report,
        })
    }

    /// One adjoint solve for dJ=sum(w_p*dp)+sum(w_Q*dQ), independent of the
    /// number of geometry, temperature and prescribed-port parameters.
    /// The symmetric grounded operator is the same one used by the primal.
    pub fn adjoint(
        &self,
        pressure_weights: &[PressureWeight],
        flow_weights: &[FlowWeight],
        cx: &Cx<'_>,
    ) -> Result<CarreauGraphGradient, CarreauGraphError> {
        checkpoint(cx)?;
        if pressure_weights.len() != self.graph.node_count()
            || flow_weights.len() != self.graph.branches.len()
        {
            return Err(GraphError::InvalidInput("Carreau adjoint weight length mismatch").into());
        }
        let mut rhs = vec![Sum::default(); self.graph.node_count()];
        for (node, weight) in pressure_weights.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            rhs[node].add(finite(weight.value(), "Carreau pressure objective weight")?)?;
        }
        for (edge, branch) in self.graph.branches.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            let value = product(
                self.conductance[edge],
                finite(flow_weights[edge].value(), "Carreau flow objective weight")?,
            )?;
            rhs[branch.from].add(value)?;
            rhs[branch.to].add(-value)?;
        }
        let (free_adjoint, linear_report) = self.solve_rhs(&rhs, cx)?;
        let mut lambda = vec![0.0; self.graph.node_count()];
        for (index, &node) in self.layout.free_nodes.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            lambda[node] = free_adjoint[index];
        }
        let edges = self.graph.branches.len();
        let mut log_radii = Vec::with_capacity(edges);
        let mut log_lengths = Vec::with_capacity(edges);
        let mut temperatures_per_k = Vec::with_capacity(edges);
        let mut boundary_gradient = vec![Sum::default(); self.graph.node_count()];
        for (node, weight) in pressure_weights.iter().enumerate() {
            periodic_checkpoint(node, cx)?;
            boundary_gradient[node].add(weight.value())?;
        }
        for (edge, branch) in self.graph.branches.iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            let corrected = finite(
                flow_weights[edge].value() - (lambda[branch.from] - lambda[branch.to]),
                "Carreau adjoint edge weight",
            )?;
            log_radii.push(product(self.responses[edge][0], corrected)?);
            log_lengths.push(product(self.responses[edge][1], corrected)?);
            temperatures_per_k.push(product(self.responses[edge][2], corrected)?);
            let value = product(self.conductance[edge], corrected)?;
            boundary_gradient[branch.from].add(value)?;
            boundary_gradient[branch.to].add(-value)?;
        }
        let mut boundary_pressures = Vec::with_capacity(self.solution.boundaries.len());
        for (index, boundary) in self.solution.boundaries.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            boundary_pressures.push(PressureWeight::new(
                boundary_gradient[boundary.node].value()?,
            ));
        }
        let mut injections = Vec::with_capacity(self.solution.injections.len());
        for (index, port) in self.solution.injections.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            injections.push(FlowWeight::new(lambda[port.node]));
        }
        checkpoint(cx)?;
        Ok(CarreauGraphGradient {
            log_radii,
            log_lengths,
            temperatures_per_k,
            boundary_pressures,
            injections,
            linear_report,
        })
    }

    fn solve_rhs(
        &self,
        node_sums: &[Sum],
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, SolveReport), CarreauGraphError> {
        let mut rhs = Vec::with_capacity(self.layout.free_nodes.len());
        for (index, &node) in self.layout.free_nodes.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            rhs.push(node_sums[node].value()?);
        }
        let operator = ConductanceLaplacian {
            graph: self.graph,
            layout: &self.layout,
            conductance: &self.conductance,
            storage_slope: None,
            inverse_diagonal: &self.inverse_diagonal,
            cx,
            failure: RefCell::new(None),
        };
        let solved = checked_cg(&operator, &operator, &rhs, self.linear, |_| {
            operator.check()
        })
        .map_err(|error| match error {
            CheckedCgError::InvalidInput(message) => {
                CarreauGraphError::Graph(GraphError::InvalidInput(message))
            }
            CheckedCgError::Interrupted(error) => error,
        })?;
        if !solved.report.converged_euclidean() {
            return Err(CarreauGraphError::LinearDidNotConverge(solved.report));
        }
        checkpoint(cx)?;
        Ok((solved.x, solved.report))
    }
}

#[cfg(test)]
mod tests;
