//! Adjoint-guided uniform exchanger sizing against a nominal outlet limit.
//!
//! Every wall is held at one declared temperature, and no external supply may
//! be colder than that wall. Under these restrictions each source-to-target
//! path contributes a nonnegative excess temperature times `exp(-s sum NTU)`;
//! fixed convex junction weights preserve monotonicity in conductance scale s.
//! Bypasses can leave a nonzero cooling floor. We therefore test feasibility at
//! the caller's bounds instead of assuming arbitrarily large exchangers work.
//!
//! This sizes frozen-flow, fixed-wall heat exchangers. It is not a geometry or
//! shared-solid optimizer, a validated engineering requirement verdict, or an
//! interval certificate of the least feasible scale.

use std::fmt;
use fs_exec::Cx;
use fs_math::det;
use fs_qty::Temperature;

use super::super::{AirSegment, BranchThermalModel, TransportError, TransportInlet,
    TransportMarch, TransportNetwork, finite, poll, positive};

/// Explicit design range, outlet requirement, and numerical work budget.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniformCoolingRequest {
    /// Hydraulic node whose mixed temperature must not exceed the limit.
    pub node: usize,
    /// Common, prescribed temperature of EVERY exchanger wall.
    pub wall_temperature: Temperature,
    /// Nominal maximum mixed temperature at the target node.
    pub outlet_limit: Temperature,
    /// Smallest admitted common multiplier of the original conductances.
    pub minimum_scale: f64,
    /// Largest admitted common multiplier of the original conductances.
    pub maximum_scale: f64,
    /// Maximum cold-side distance from the limit at an interior root, K.
    pub temperature_tolerance_k: f64,
    /// Maximum natural-log width of the failing/passing sizing bracket.
    pub log_scale_tolerance: f64,
    /// Maximum primal marches, including both endpoint evaluations.
    pub max_evaluations: usize,
}

/// A nominally passing design and its last numerical search bracket.
#[derive(Debug, Clone, PartialEq)]
pub struct UniformCoolingDesign {
    /// Common multiplier to apply with `scaled_conductance`.
    pub scale: f64,
    /// Target mixed temperature from the returned design's actual march.
    pub outlet_temperature: Temperature,
    /// Lower bracket scale; nominally failing unless `at_lower_bound` is true.
    pub lower_scale: f64,
    /// Target temperature at the retained lower scale.
    pub lower_outlet_temperature: Temperature,
    /// Derivative of target temperature with respect to log(common scale), K.
    pub slope_per_log_scale_k: f64,
    /// Primal march count; an adjoint sweep does not consume a primal evaluation.
    pub evaluations: usize,
    /// True when the smallest declared scale already meets the limit.
    pub at_lower_bound: bool,
    /// The full admitted transport result for the passing design, not a proxy.
    pub march: TransportMarch,
}

/// Typed configuration, physical-admission and bounded-search failures.
#[derive(Debug)]
pub enum CoolingSizingError {
    /// Underlying hydraulic/thermal/cancellation failure, preserved intact.
    Transport(TransportError),
    /// Unusable request or a numerically unresolved bracket.
    InvalidInput(&'static str),
    /// A supply is colder than the common wall; monotone cooling is not assured.
    SupplyBelowWall { node: usize },
    /// The requested node does not have a transported temperature.
    UnknownTarget { node: usize },
    /// Even the largest admitted design does not satisfy the limit.
    Unattainable { maximum_scale: f64, outlet_temperature: Temperature },
    /// The work budget ended before both requested numerical tolerances held.
    BudgetExhausted { evaluations: usize, lower_scale: f64, upper_scale: f64 },
}
impl From<TransportError> for CoolingSizingError {
    fn from(error: TransportError) -> Self { Self::Transport(error) }
}
impl fmt::Display for CoolingSizingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "uniform cooling sizing: {self:?}") }
}
impl std::error::Error for CoolingSizingError {}

struct Sample { scale: f64, log_scale: f64, temperature: f64, slope: f64, march: TransportMarch }

impl<'flow> TransportNetwork<'flow> {
    /// Rebuild the same admitted network with every h multiplied by `scale`.
    /// Areas, flows, temperatures, property values and transport budgets remain
    /// unchanged. Stored actual-flow row order is converted back to declared
    /// edge order before admission, so reverse-flow segments are not reversed
    /// twice. This also supplies the candidate-construction seam for optimizers.
    pub fn scaled_conductance(&self, cx: &Cx<'_>, scale: f64) -> Result<Self, TransportError> {
        poll(cx)?;
        positive(scale, "conductance scale must be finite and positive")?;
        let mut models = Vec::with_capacity(self.models.len());
        for (branch, model) in self.models.iter().enumerate() {
            poll(cx)?;
            models.push(match model {
                BranchThermalModel::Adiabatic => BranchThermalModel::Adiabatic,
                BranchThermalModel::Exchange(segments) => {
                    let mut rows = Vec::with_capacity(segments.len());
                    for segment in segments {
                        poll(cx)?;
                        rows.push(AirSegment::new(segment.region(), segment.area_m2(),
                            segment.htc_w_per_m2_k() * scale)?);
                    }
                    if self.flow.branches[branch].flow.value() < 0.0 { rows.reverse(); }
                    BranchThermalModel::Exchange(rows)
                }
            });
        }
        let mut inlets = Vec::new();
        for (node, temperature) in self.inlet_temperatures.iter().enumerate() {
            poll(cx)?;
            if let Some(value) = temperature { inlets.push(TransportInlet { node, temperature: Temperature::new(*value) }); }
        }
        let network = Self::new(cx, self.flow, self.air, models, &inlets, self.config)?;
        match &self.feedback {
            Some(feedback) => network.with_recirculation(cx, feedback.links.clone(), feedback.tolerance_k),
            None => Ok(network),
        }
    }

    /// Find the least nominally passing uniform conductance scale within the
    /// declared range, to the requested numerical bracket width. Safeguarded
    /// Newton uses one adjoint per candidate; bisection preserves the bracket.
    /// Success always returns the PASSING side, never an infeasible midpoint.
    /// The lower-bound case makes no claim of feasibility below that bound.
    pub fn size_uniform_cooling(&self, cx: &Cx<'_>, request: UniformCoolingRequest)
        -> Result<UniformCoolingDesign, CoolingSizingError>
    {
        poll(cx)?;
        for (name, value) in [("sizing wall temperature", request.wall_temperature.value()),
            ("sizing outlet limit", request.outlet_limit.value()), ("minimum conductance scale", request.minimum_scale),
            ("maximum conductance scale", request.maximum_scale), ("sizing kelvin tolerance", request.temperature_tolerance_k),
            ("sizing logarithmic tolerance", request.log_scale_tolerance)] { positive(value, name)?; }
        if request.minimum_scale >= request.maximum_scale || request.max_evaluations < 2 {
            return Err(CoolingSizingError::InvalidInput("ordered distinct scale bounds and at least two evaluations are required"));
        }
        if request.node >= self.order.len() { return Err(CoolingSizingError::UnknownTarget { node: request.node }); }
        for (node, temperature) in self.inlet_temperatures.iter().enumerate() {
            poll(cx)?;
            if temperature.is_some_and(|t| t < request.wall_temperature.value()) {
                return Err(CoolingSizingError::SupplyBelowWall { node });
            }
        }
        let walls = vec![request.wall_temperature.value(); self.offsets.last().copied().unwrap_or(0)];
        let evaluate = |scale: f64| -> Result<Sample, CoolingSizingError> {
            let network = self.scaled_conductance(cx, scale)?;
            let linearization = network.linearize(cx, &walls)?;
            let temperature = linearization.primal.node_temperatures_k[request.node]
                .ok_or(CoolingSizingError::UnknownTarget { node: request.node })?;
            let mut objective = linearization.zero_objective(); objective.node_temperatures[request.node] = 1.0;
            let gradient = linearization.pullback(cx, &objective)?;
            let mut slope = 0.0;
            for value in gradient.log_conductances { poll(cx)?; slope = finite(slope + value, "uniform sizing derivative")?; }
            Ok(Sample { scale, log_scale: finite(det::ln(scale), "logarithmic conductance scale")?,
                temperature, slope, march: linearization.primal })
        };
        let limit = request.outlet_limit.value();
        let mut lower = evaluate(request.minimum_scale)?;
        if lower.temperature <= limit {
            poll(cx)?;
            return Ok(design(lower.scale, lower.temperature, lower, 1, true));
        }
        let mut upper = evaluate(request.maximum_scale)?;
        if upper.temperature > limit {
            return Err(CoolingSizingError::Unattainable { maximum_scale: upper.scale,
                outlet_temperature: Temperature::new(upper.temperature) });
        }
        let mut evaluations = 2;
        loop {
            poll(cx)?;
            let width = upper.log_scale - lower.log_scale;
            if width <= request.log_scale_tolerance && limit - upper.temperature <= request.temperature_tolerance_k {
                return Ok(design(lower.scale, lower.temperature, upper, evaluations, false));
            }
            if evaluations == request.max_evaluations {
                return Err(CoolingSizingError::BudgetExhausted { evaluations, lower_scale: lower.scale, upper_scale: upper.scale });
            }
            let proposed = if upper.slope < 0.0 {
                upper.log_scale - (upper.temperature - limit) / upper.slope
            } else { f64::NAN };
            let next_log = if proposed.is_finite() && proposed > lower.log_scale + 0.1*width
                && proposed < upper.log_scale - 0.1*width { proposed }
                else { 0.5*lower.log_scale + 0.5*upper.log_scale };
            let scale = det::exp(next_log);
            if !(scale.is_finite() && scale > lower.scale && scale < upper.scale) {
                return Err(CoolingSizingError::InvalidInput("requested sizing precision exceeds representable bracket resolution"));
            }
            let next = evaluate(scale)?; evaluations += 1;
            if next.temperature <= limit { upper = next; } else { lower = next; }
        }
    }
}

fn design(lower_scale: f64, lower_temperature: f64, sample: Sample, evaluations: usize, at_lower_bound: bool)
    -> UniformCoolingDesign
{
    UniformCoolingDesign { scale: sample.scale, outlet_temperature: Temperature::new(sample.temperature),
        lower_scale, lower_outlet_temperature: Temperature::new(lower_temperature),
        slope_per_log_scale_k: sample.slope, evaluations, at_lower_bound, march: sample.march }
}

#[cfg(test)]
mod tests;
