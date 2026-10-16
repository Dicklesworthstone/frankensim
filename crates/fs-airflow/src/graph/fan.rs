//! Fan operating points for a fixed passive graph. These are estimates, not
//! the interval-certified series/parallel operating points in the parent crate.

use core::fmt;

use fs_exec::Cx;
use fs_qty::{Pressure, VolumetricFlowRate};

use super::{FixedPressure, FlowState, GraphError, GraphSolution, GraphSolveConfig, LossGraph, checkpoint, finite};
use crate::{AirflowError, FanBank};

/// Explicit scalar operating-point budgets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GraphFanConfig {
    /// Maximum bisection iterations in the retained fan-flow domain.
    pub max_iterations: usize,
    /// Maximum absolute flow-bracket width and final reconstruction error.
    pub absolute_flow_tolerance: VolumetricFlowRate,
    /// Maximum absolute fan/system pressure mismatch, in Pa.
    pub absolute_pressure_tolerance: Pressure,
}

impl GraphFanConfig {
    fn validate(self) -> Result<(), GraphFanError> {
        let flow = self.absolute_flow_tolerance.value();
        let pressure = self.absolute_pressure_tolerance.value();
        if self.max_iterations == 0 || !(flow.is_finite() && flow > 0.0)
            || !(pressure.is_finite() && pressure > 0.0)
        {
            return Err(GraphFanError::InvalidInput("fan iteration count and finite absolute tolerances must be positive"));
        }
        Ok(())
    }
}

/// Nominal fan/graph intersection with independently recomputed residuals.
#[derive(Debug, Clone, PartialEq)]
pub struct GraphOperatingPoint {
    /// Retained fan curve, uncertainty, source, arrangement, and speed.
    pub fan: FanBank,
    /// Pressure and flow solution at the operating point.
    pub network: GraphSolution,
    /// Recomputed flow supplied at the graph inlet, in m^3/s.
    pub flow: VolumetricFlowRate,
    /// Fan pressure evaluated at the recomputed flow, in Pa.
    pub fan_pressure: Pressure,
    /// Actual graph inlet-to-outlet pressure, in Pa.
    pub network_pressure: Pressure,
    /// Signed fan pressure minus network pressure, in Pa.
    pub pressure_residual: Pressure,
    /// Completed scalar bisections; graph sweeps are recorded separately.
    pub iterations: usize,
}

/// Fan, graph, retained-domain, or convergence refusal.
#[derive(Debug)]
pub enum GraphFanError {
    /// Invalid terminals, scaling, or scalar configuration.
    InvalidInput(&'static str),
    /// A graph solve, arithmetic check, or context checkpoint refused.
    Graph(GraphError),
    /// The existing fan model refused evaluation.
    Fan(AirflowError),
    /// There is no connected through-path between the selected fan terminals.
    NoThroughPath,
    /// No nominal intersection lies in the retained admissible fan domain.
    OutsideFanDomain(&'static str),
    /// The bounded scalar solve or its reconstructed result missed tolerance.
    DidNotConverge {
        /// Completed scalar iterations.
        iterations: usize,
        /// Last pressure mismatch, in Pa.
        pressure_residual: Pressure,
    },
}

impl From<GraphError> for GraphFanError {
    fn from(error: GraphError) -> Self { Self::Graph(error) }
}

impl From<AirflowError> for GraphFanError {
    fn from(error: AirflowError) -> Self { Self::Fan(error) }
}

impl fmt::Display for GraphFanError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "airflow graph/fan: {self:?}")
    }
}

impl std::error::Error for GraphFanError {}

impl LossGraph {
    /// Intersect an arbitrary passive graph with an existing identical-fan bank.
    ///
    /// The pressure outlet is the zero reference. Fan-law scaling, arrangement,
    /// and stall bounds come from `FanBank`; no extrapolation is performed.
    ///
    /// A fixed quadratic graph is homogeneous: scaling every nodal pressure by
    /// `s^2` scales every flow by `s`. The reference graph is solved at the
    /// largest retained fan pressure, not at an arbitrary 1 Pa. Subsequent
    /// pressure scaling is at most one, so it cannot amplify the reference
    /// graph's absolute conservation tolerance. Final flows, conservation, and
    /// fan mismatch are nevertheless recomputed and checked before returning.
    ///
    /// This numerical reduction is NOT a certified equivalent resistance or
    /// uncertainty enclosure. It cannot be used to mint `OperatingPoint`'s
    /// interval-Newton certificate. All graph and fan metadata stays retained.
    pub fn solve_with_fan(
        &self, inlet: usize, outlet: usize, fan: &FanBank,
        graph_config: GraphSolveConfig, fan_config: GraphFanConfig, cx: &Cx<'_>,
    ) -> Result<GraphOperatingPoint, GraphFanError> {
        checkpoint(cx)?;
        graph_config.validate()?;
        fan_config.validate()?;
        self.check_fan_terminals(inlet, outlet, cx)?;
        let factor = fan.flow_factor();
        let pressure_factor = fan.pressure_factor();
        if !(factor.is_finite() && factor > 0.0 && pressure_factor.is_finite() && pressure_factor > 0.0) {
            return Err(GraphFanError::InvalidInput("fan-law scaling is not representable"));
        }
        let minimum = finite(fan.curve().admissible_min_flow().value() * factor, "fan minimum flow")?;
        let maximum = finite(fan.curve().points().last().expect("validated fan curve").flow.value() * factor, "fan maximum flow")?;
        if minimum >= maximum {
            return Err(GraphFanError::InvalidInput("scaled fan domain has no representable width"));
        }
        let peak = finite(fan.pressure_at(VolumetricFlowRate::new(minimum))?.value(), "peak fan pressure")?;
        if peak < 0.0 {
            return Err(GraphFanError::InvalidInput("fan pressure must be nonnegative"));
        }
        if peak == 0.0 {
            if minimum != 0.0 {
                return Err(GraphFanError::OutsideFanDomain("zero pressure cannot supply the minimum admitted flow"));
            }
            let network = self.solve(&[fixed(inlet, 0.0), fixed(outlet, 0.0)], graph_config, cx)?;
            return checked_operating_point(network, inlet, outlet, fan, 0.0, 0, fan_config, cx);
        }
        let reference = self.solve(&[fixed(inlet, peak), fixed(outlet, 0.0)], graph_config, cx)?;
        let reference_flow = reference.node_outflows[inlet].value();
        if !(reference_flow.is_finite() && reference_flow > 0.0) {
            return Err(GraphFanError::InvalidInput("reference graph through-flow is not representable"));
        }
        if reference_flow < minimum {
            return Err(GraphFanError::OutsideFanDomain("intersection is below the declared minimum fan flow"));
        }
        let mut low = minimum;
        let mut high = maximum.min(reference_flow);
        let residual = |flow: f64| -> Result<f64, GraphFanError> {
            let ratio = flow / reference_flow;
            let system_pressure = finite(peak * ratio * ratio, "system pressure")?;
            Ok(finite(fan.pressure_at(VolumetricFlowRate::new(flow))?.value() - system_pressure, "fan/system mismatch")?)
        };
        let low_residual = residual(low)?;
        let high_residual = residual(high)?;
        if low_residual < 0.0 || high_residual > 0.0 {
            return Err(GraphFanError::OutsideFanDomain("no sign bracket inside the retained fan curve"));
        }
        let (flow, iterations) = if low_residual == 0.0 {
            (low, 0)
        } else if high_residual == 0.0 {
            (high, 0)
        } else {
            let mut found = None;
            let mut last_residual = low_residual;
            let mut completed = 0;
            for iteration in 1..=fan_config.max_iterations {
                checkpoint(cx)?;
                completed = iteration;
                let middle = 0.5 * low + 0.5 * high;
                let value = residual(middle)?;
                last_residual = value;
                if value == 0.0 || (high - low <= fan_config.absolute_flow_tolerance.value()
                    && value.abs() <= fan_config.absolute_pressure_tolerance.value())
                {
                    found = Some((middle, iteration));
                    break;
                }
                if middle <= low || middle >= high { break; }
                if value > 0.0 { low = middle; } else { high = middle; }
            }
            found.ok_or(GraphFanError::DidNotConverge {
                iterations: completed,
                pressure_residual: Pressure::new(last_residual),
            })?
        };
        let ratio = flow / reference_flow;
        let pressures: Vec<f64> = reference.pressures.iter()
            .map(|pressure| pressure.value() * ratio * ratio).collect();
        let boundaries = [fixed(inlet, pressures[inlet]), fixed(outlet, pressures[outlet])];
        let prescribed = self.fixed_pressures(&boundaries, cx)?;
        let mut state = FlowState::new(self.node_count(), self.branches.len());
        if !self.evaluate(&pressures, &prescribed, graph_config, cx, &mut state)? {
            return Err(state.nonconvergence(reference.sweeps).into());
        }
        let network = self.finish(pressures, &boundaries, state, reference.sweeps, cx)?;
        checked_operating_point(network, inlet, outlet, fan, flow, iterations, fan_config, cx)
    }

    fn check_fan_terminals(&self, inlet: usize, outlet: usize, cx: &Cx<'_>) -> Result<(), GraphFanError> {
        if inlet >= self.node_count() || outlet >= self.node_count() || inlet == outlet {
            return Err(GraphFanError::InvalidInput("fan terminals must be distinct valid nodes"));
        }
        let mut seen = vec![false; self.node_count()];
        let mut stack = vec![inlet];
        seen[inlet] = true;
        while let Some(node) = stack.pop() {
            checkpoint(cx)?;
            if node == outlet { return Ok(()); }
            for (ordinal, &edge) in self.adjacency[node].iter().enumerate() {
                super::periodic_checkpoint(ordinal, cx)?;
                let other = self.other(edge, node);
                if !seen[other] { seen[other] = true; stack.push(other); }
            }
        }
        Err(GraphFanError::NoThroughPath)
    }
}

fn fixed(node: usize, pressure: f64) -> FixedPressure {
    FixedPressure { node, pressure: Pressure::new(pressure) }
}

#[allow(clippy::too_many_arguments)]
fn checked_operating_point(
    network: GraphSolution, inlet: usize, outlet: usize, fan: &FanBank, trial_flow: f64,
    iterations: usize, config: GraphFanConfig, cx: &Cx<'_>,
) -> Result<GraphOperatingPoint, GraphFanError> {
    let flow = network.node_outflows[inlet];
    let fan_pressure = fan.pressure_at(flow)?;
    let network_pressure = Pressure::new(finite(
        network.pressures[inlet].value() - network.pressures[outlet].value(), "final graph pressure",
    )?);
    let difference = finite(fan_pressure.value() - network_pressure.value(), "final fan mismatch")?;
    if difference.abs() > config.absolute_pressure_tolerance.value()
        || (flow.value() - trial_flow).abs() > config.absolute_flow_tolerance.value()
    {
        return Err(GraphFanError::DidNotConverge { iterations, pressure_residual: Pressure::new(difference) });
    }
    checkpoint(cx)?;
    Ok(GraphOperatingPoint {
        fan: fan.clone(), network, flow, fan_pressure, network_pressure,
        pressure_residual: Pressure::new(difference), iterations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{FanArrangement, FanCurve, FanPoint, SourceProvenance, ToleranceBasis};
    use super::super::tests::{bridge, close, config, edge, with_cx};

    fn fan(count: usize, arrangement: FanArrangement, speed: f64, stall: f64, end: f64) -> FanBank {
        FanBank::new(FanCurve::new(
            "linear fixture",
            vec![FanPoint::new(VolumetricFlowRate::new(0.0), Pressure::new(100.0)),
                FanPoint::new(VolumetricFlowRate::new(end), Pressure::new(100.0 - 10.0 * end))],
            SourceProvenance::new("analytic curve", "linear-fan-v1"), 0.1,
            ToleranceBasis::EngineeringAllowance, VolumetricFlowRate::new(stall), (0.25, 2.0),
        ).unwrap(), count, arrangement, speed).unwrap()
    }

    fn fan_config() -> GraphFanConfig {
        GraphFanConfig { max_iterations: 128, absolute_flow_tolerance: VolumetricFlowRate::new(1.0e-8),
            absolute_pressure_tolerance: Pressure::new(1.0e-7) }
    }

    #[test]
    fn g1_fan_arrangement_and_speed_match_quadratic_closed_forms() {
        let graph = LossGraph::new(2, vec![edge("loss", 0, 1, 4.0)]).unwrap();
        for (count, arrangement, speed, intercept, slope) in [
            (1, FanArrangement::Series, 1.0, 100.0, 10.0),
            (2, FanArrangement::Series, 1.0, 200.0, 20.0),
            (2, FanArrangement::Parallel, 1.0, 100.0, 5.0),
            (1, FanArrangement::Series, 0.5, 25.0, 5.0),
        ] {
            let bank = fan(count, arrangement, speed, 0.0, 10.0);
            let expected = (-slope + (slope * slope + 16.0_f64 * intercept).sqrt()) / 8.0;
            let point = with_cx(|cx| graph.solve_with_fan(0, 1, &bank, config(), fan_config(), cx)).unwrap();
            close(point.flow.value(), expected);
            close(point.network_pressure.value(), 4.0 * expected * expected);
            assert!(point.pressure_residual.value().abs() <= fan_config().absolute_pressure_tolerance.value());
            assert_eq!(point.fan, bank);
        }
    }

    #[test]
    fn g1_fan_drives_the_non_series_parallel_bridge() {
        let graph = bridge();
        let bank = fan(1, FanArrangement::Series, 1.0, 0.0, 10.0);
        let point = with_cx(|cx| graph.solve_with_fan(0, 3, &bank, config(), fan_config(), cx)).unwrap();
        let resistance = 25.0 / 16.0;
        let expected = (-10.0 + (100.0_f64 + 400.0 * resistance).sqrt()) / (2.0 * resistance);
        close(point.flow.value(), expected);
        close(point.network.branches[2].flow.value(), -expected / 4.0);
    }

    #[test]
    fn g0_stall_and_retained_curve_limits_are_not_extrapolated() {
        let graph = LossGraph::new(2, vec![edge("loss", 0, 1, 4.0)]).unwrap();
        for bank in [fan(1, FanArrangement::Series, 1.0, 4.0, 10.0),
            fan(1, FanArrangement::Series, 1.0, 0.0, 2.0)] {
            let result = with_cx(|cx| graph.solve_with_fan(0, 1, &bank, config(), fan_config(), cx));
            assert!(matches!(result, Err(GraphFanError::OutsideFanDomain(_))));
        }
    }

    #[test]
    fn g0_fan_requires_a_through_path() {
        let graph = LossGraph::new(4, vec![edge("a", 0, 1, 1.0), edge("b", 2, 3, 1.0)]).unwrap();
        let bank = fan(1, FanArrangement::Series, 1.0, 0.0, 10.0);
        let result = with_cx(|cx| graph.solve_with_fan(0, 3, &bank, config(), fan_config(), cx));
        assert!(matches!(result, Err(GraphFanError::NoThroughPath)));
    }

    #[test]
    fn g0_scalar_exhaustion_does_not_return_an_operating_point() {
        let graph = LossGraph::new(2, vec![edge("loss", 0, 1, 4.0)]).unwrap();
        let bank = fan(1, FanArrangement::Series, 1.0, 0.0, 10.0);
        let mut budget = fan_config();
        budget.max_iterations = 1;
        let result = with_cx(|cx| graph.solve_with_fan(0, 1, &bank, config(), budget, cx));
        assert!(matches!(result, Err(GraphFanError::DidNotConverge { .. })));
    }
}
