//! Local derivatives of a fan-driven graph, including operating-point feedback.
//!
//! The inlet pressure is NOT held fixed. For `p = F(Q)` the scalar feedback
//! denominator is `1 - F_Q Q_p`; its network derivatives come from the same
//! accepted hydraulic state as the objective. Fan-law speed and retained
//! pressure samples are controls, while fan count and arrangement stay fixed.
//! Knots, retained-domain boundaries, and stagnant branches refuse. These are
//! nominal within-segment derivatives, not certificates or stall sensitivities.

use core::fmt;

use fs_exec::Cx;
use fs_qty::{Pressure, VolumetricFlowRate};
use fs_solver::SolveReport;

use super::{
    FlowWeight, GraphGradient, GraphLinearization, GraphSensitivityConfig,
    GraphSensitivityError, GraphTangent, PressureWeight,
};
use super::super::{
    GraphError, GraphFanConfig, GraphFanError, GraphSolution, GraphSolveConfig,
    LossGraph, checkpoint, finite, periodic_checkpoint,
};
use crate::FanBank;

/// An operating-point, hydraulic, or local fan-derivative refusal.
#[derive(Debug)]
pub enum GraphFanSensitivityError {
    /// The existing fan/network operating-point solver refused.
    OperatingPoint(GraphFanError),
    /// The hydraulic derivative model or solve refused.
    Hydraulic(GraphSensitivityError),
    /// A two-sided local fan derivative is not admitted here.
    NonDifferentiableFan(&'static str),
    /// Recomputed differentiated fan balance missed the linear tolerance.
    FeedbackDidNotConverge {
        /// Signed `dp - F_Q dQ - dF_control`, in Pa.
        pressure_residual: Pressure,
        /// Residual divided by the largest absolute balance term.
        relative_residual: f64,
    },
}

impl From<GraphFanError> for GraphFanSensitivityError {
    fn from(error: GraphFanError) -> Self { Self::OperatingPoint(error) }
}
impl From<GraphSensitivityError> for GraphFanSensitivityError {
    fn from(error: GraphSensitivityError) -> Self { Self::Hydraulic(error) }
}
impl From<GraphError> for GraphFanSensitivityError {
    fn from(error: GraphError) -> Self { Self::Hydraulic(error.into()) }
}
impl fmt::Display for GraphFanSensitivityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "airflow fan sensitivity: {self:?}")
    }
}
impl std::error::Error for GraphFanSensitivityError {}

/// Total directional derivative, with the fan/network balance closed.
#[derive(Debug, Clone)]
pub struct GraphFanTangent {
    /// Pressure and branch-flow derivatives, including reservoir exchanges.
    pub network: GraphTangent,
    /// Total fan-flow derivative (m^3/s).
    pub fan_flow: VolumetricFlowRate,
    /// Fan/network pressure derivative (Pa); outlet reference stays zero.
    pub fan_pressure: Pressure,
    /// Recomputed differentiated fan-balance residual, not an error enclosure.
    pub pressure_residual: Pressure,
}

/// Gradient of a dimensionless objective at the moving operating point.
#[derive(Debug, Clone)]
pub struct GraphFanGradient {
    /// Common resistance-scale derivative in edge order, including fan feedback.
    /// A mixed edge scales both its linear and quadratic coefficients together.
    pub log_resistances: Vec<f64>,
    /// `dJ/d(log speed_ratio)`; divide by the speed ratio for `dJ/dspeed_ratio`.
    pub log_speed_ratio: f64,
    /// `dJ/d(base-curve pressure)` in retained sample order (1/Pa).
    pub curve_pressures: Vec<PressureWeight>,
    /// Objective's hydraulic adjoint acceptance report.
    pub objective_linear_report: SolveReport,
    /// Reused inlet-flow adjoint acceptance report.
    pub inlet_linear_report: SolveReport,
}

#[derive(Debug)]
struct FanSegment {
    index: usize,
    fraction: f64,
    // Pa per (m^3/s) of total bank flow, at the current speed.
    slope: f64,
}

/// Immutable accepted fan/network state, reusable across objectives.
#[derive(Debug)]
pub struct GraphFanLinearization<'g> {
    hydraulic: GraphLinearization<'g>,
    fan: FanBank,
    inlet: usize,
    segment: FanSegment,
    // Partial fan pressure derivative w.r.t. log speed at fixed total flow.
    speed_pressure: f64,
    // Gradient of Q/(1 m^3/s); numeric values rescale to SI flow derivatives.
    inlet_gradient: GraphGradient,
    feedback_denominator: f64,
}

impl LossGraph {
    /// Solve a fan-driven graph and admit total operating-point derivatives.
    ///
    /// Reuses the accepted operating point without resolving at a frozen inlet
    /// pressure. A full fan-flow tolerance (and floating-point roundoff margin)
    /// must separate the nominal flow from every knot and domain endpoint.
    /// Speed must be strictly inside its declared domain. These are local
    /// admission rules, not a certified distance to a true root or a kink.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_fan_linearized(
        &self,
        inlet: usize,
        outlet: usize,
        fan: &FanBank,
        graph_config: GraphSolveConfig,
        fan_config: GraphFanConfig,
        sensitivity: GraphSensitivityConfig,
        cx: &Cx<'_>,
    ) -> Result<GraphFanLinearization<'_>, GraphFanSensitivityError> {
        checkpoint(cx)?;
        sensitivity.validate()?;
        let point = self.solve_with_fan(inlet, outlet, fan, graph_config, fan_config, cx)?;
        let segment = admit_segment(fan, point.flow.value(), fan_config, cx)?;
        let speed_pressure = finite(
            2.0 * point.fan_pressure.value() - point.flow.value() * segment.slope,
            "fan log-speed partial",
        )?;
        let hydraulic = self.linearize_accepted(point.network, sensitivity, cx)?;
        let mut weights = Vec::with_capacity(self.branches().len());
        for (edge, branch) in self.branches().iter().enumerate() {
            periodic_checkpoint(edge, cx)?;
            let value = if branch.from == inlet { 1.0 }
                else if branch.to == inlet { -1.0 } else { 0.0 };
            weights.push(FlowWeight::new(value));
        }
        let inlet_gradient = hydraulic.adjoint(
            &vec![PressureWeight::new(0.0); self.node_count()], &weights, cx,
        )?;
        // solve_with_fan owns the boundary order: inlet, then zero-reference outlet.
        let flow_per_pressure = inlet_gradient.boundary_pressures[0].value();
        let feedback_denominator = finite(
            1.0 - segment.slope * flow_per_pressure, "fan feedback denominator",
        )?;
        if flow_per_pressure <= 0.0 || feedback_denominator < 1.0 {
            return Err(GraphFanSensitivityError::NonDifferentiableFan(
                "positive inlet conductance and a nonincreasing fan curve required",
            ));
        }
        checkpoint(cx)?;
        Ok(GraphFanLinearization {
            hydraulic, fan: point.fan, inlet, segment, speed_pressure,
            inlet_gradient, feedback_denominator,
        })
    }
}

impl GraphFanLinearization<'_> {
    /// Accepted nominal network, retaining its physical metadata unchanged.
    #[must_use]
    pub fn solution(&self) -> &GraphSolution { self.hydraulic.solution() }

    /// Retained fan curve, arrangement and speed.
    #[must_use]
    pub fn fan(&self) -> &FanBank { &self.fan }

    /// Total derivative for resistance, log-speed and base-curve pressure controls.
    ///
    /// Resistance directions are relative to nominal R. Curve directions are
    /// pressures in original sample order; knot flow coordinates stay fixed.
    pub fn tangent(
        &self,
        log_resistance_direction: &[f64],
        log_speed_ratio_direction: f64,
        curve_pressure_direction: &[Pressure],
        cx: &Cx<'_>,
    ) -> Result<GraphFanTangent, GraphFanSensitivityError> {
        checkpoint(cx)?;
        if curve_pressure_direction.len() != self.fan.curve().points().len() {
            return Err(GraphSensitivityError::InvalidInput("fan pressure direction length mismatch").into());
        }
        finite(log_speed_ratio_direction, "fan log-speed direction")?;
        for (index, value) in curve_pressure_direction.iter().enumerate() {
            periodic_checkpoint(index, cx)?;
            finite(value.value(), "fan pressure direction")?;
        }
        let k = self.segment.index;
        let t = self.segment.fraction;
        let control = finite(self.speed_pressure * log_speed_ratio_direction
            + self.fan.pressure_factor() * ((1.0 - t) * curve_pressure_direction[k].value()
                + t * curve_pressure_direction[k + 1].value()), "fan control pressure")?;
        let fixed = self.hydraulic.tangent(log_resistance_direction, &[Pressure::new(0.0); 2], cx)?;
        let pressure = finite((self.segment.slope * fixed.node_outflows[self.inlet].value()
            + control) / self.feedback_denominator, "moving inlet pressure")?;
        // Re-solve the complete tangent RHS, so its report belongs to the final
        // network derivative rather than to a superposition of two estimates.
        let network = self.hydraulic.tangent(
            log_resistance_direction, &[Pressure::new(pressure), Pressure::new(0.0)], cx,
        )?;
        let flow = network.node_outflows[self.inlet];
        let feedback = finite(self.segment.slope * flow.value(), "fan flow feedback")?;
        let residual = finite(pressure - feedback - control, "differentiated fan balance")?;
        let scale = pressure.abs().max(feedback.abs()).max(control.abs());
        let relative = if scale == 0.0 { 0.0 } else { residual.abs() / scale };
        if relative > self.hydraulic.config.linear.tolerance {
            return Err(GraphFanSensitivityError::FeedbackDidNotConverge {
                pressure_residual: Pressure::new(residual), relative_residual: relative,
            });
        }
        checkpoint(cx)?;
        Ok(GraphFanTangent {
            network, fan_flow: flow, fan_pressure: Pressure::new(pressure),
            pressure_residual: Pressure::new(residual),
        })
    }

    /// Total adjoint gradient, including the movement of the inlet pressure.
    ///
    /// Caller weights and explicit objective dependence follow `GraphLinearization::adjoint`.
    /// Only the active segment's pressure samples have nonzero curve gradients.
    pub fn adjoint(
        &self, pressure_weights: &[PressureWeight], flow_weights: &[FlowWeight], cx: &Cx<'_>,
    ) -> Result<GraphFanGradient, GraphFanSensitivityError> {
        checkpoint(cx)?;
        let objective = self.hydraulic.adjoint(pressure_weights, flow_weights, cx)?;
        let multiplier = finite(objective.boundary_pressures[0].value()
            / self.feedback_denominator, "fan adjoint multiplier")?;
        let influence = finite(multiplier * self.segment.slope, "fan resistance feedback")?;
        let mut log_resistances = Vec::with_capacity(objective.log_resistances.len());
        for (edge, (&direct, &flow)) in objective.log_resistances.iter()
            .zip(&self.inlet_gradient.log_resistances).enumerate()
        {
            periodic_checkpoint(edge, cx)?;
            log_resistances.push(finite(direct + influence * flow, "total resistance gradient")?);
        }
        let log_speed_ratio = finite(multiplier * self.speed_pressure, "fan speed gradient")?;
        let mut curve_pressures = vec![PressureWeight::new(0.0); self.fan.curve().points().len()];
        for (offset, weight) in [1.0 - self.segment.fraction, self.segment.fraction].into_iter().enumerate() {
            curve_pressures[self.segment.index + offset] = PressureWeight::new(finite(
                multiplier * self.fan.pressure_factor() * weight, "fan curve pressure gradient",
            )?);
        }
        checkpoint(cx)?;
        Ok(GraphFanGradient {
            log_resistances, log_speed_ratio, curve_pressures,
            objective_linear_report: objective.linear_report,
            inlet_linear_report: self.inlet_gradient.linear_report.clone(),
        })
    }
}

fn admit_segment(
    fan: &FanBank, flow: f64, config: GraphFanConfig, cx: &Cx<'_>,
) -> Result<FanSegment, GraphFanSensitivityError> {
    let (low_speed, high_speed) = fan.curve().speed_ratio_domain;
    if !(fan.speed_ratio() > low_speed && fan.speed_ratio() < high_speed) {
        return Err(GraphFanSensitivityError::NonDifferentiableFan("speed lies on a retained-domain boundary"));
    }
    let factor = fan.flow_factor();
    let margin = config.absolute_flow_tolerance.value().max(8.0 * f64::EPSILON * flow.abs());
    let minimum = finite(fan.curve().admissible_min_flow().value() * factor, "fan derivative domain")?;
    if flow - minimum <= margin {
        return Err(GraphFanSensitivityError::NonDifferentiableFan("flow is too close to the stall/domain boundary"));
    }
    for (index, pair) in fan.curve().points().windows(2).enumerate() {
        periodic_checkpoint(index, cx)?;
        let low = finite(pair[0].flow.value() * factor, "fan segment flow")?;
        let high = finite(pair[1].flow.value() * factor, "fan segment flow")?;
        if flow >= low && flow <= high {
            if flow - low <= margin || high - flow <= margin {
                return Err(GraphFanSensitivityError::NonDifferentiableFan("flow is too close to a curve knot/domain endpoint"));
            }
            let difference = pair[1].pressure.value() - pair[0].pressure.value();
            let slope = finite(difference * fan.pressure_factor() / (high - low), "fan segment slope")?;
            if slope > 0.0 || (difference != 0.0 && slope == 0.0) {
                return Err(GraphFanSensitivityError::NonDifferentiableFan("fan slope is not representable"));
            }
            return Ok(FanSegment { index, fraction: (flow - low) / (high - low), slope });
        }
    }
    Err(GraphFanSensitivityError::NonDifferentiableFan("no interior retained fan segment"))
}

#[cfg(test)]
mod tests;
