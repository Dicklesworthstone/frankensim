//! Source-control derivatives through the actual implicit scan endpoint.
//!
//! The existing enthalpy/radiation owner revalidates the primal and builds the
//! complete physical tangent. One transpose solve supplies every source's
//! center/width/watts pullback plus the previous enthalpy seed for a reverse
//! sweep on a FIXED time grid. No solver iterations, adaptive decisions, event
//! times, material selections or source-support admission rules are differentiated.

/// Full fixed-grid inverse objectives with authored-schedule gradients.
pub mod objective;

use super::{ScanningEnthalpyRun, ScanStep, events, natural_only, patches_before, poll};
use super::super::super::moving::sensitivity::{MovingHeatGradient, MovingHeatLinearization};
use crate::{ConductionError, ConductionProblem, LinearConfig};
use crate::transient::enthalpy::EnthalpyStepConfig;
use crate::transient::enthalpy::adjoint::EnthalpyAdjointError;
use crate::transient::enthalpy::heterogeneous::{HeterogeneousEnthalpyError,
    HeterogeneousEnthalpyStepGradient, HeterogeneousEnthalpyStepLinearization};
use fs_exec::Cx;

/// Extra retained-source and radiation-tangent caps for one linearized step.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanSensitivityBudget {
    /// Maximum sources times mesh vertices. Each binding retains at most three
    /// nodal scalar arrays; this count is not a byte/allocator-overhead bound.
    pub max_source_vertices: usize,
    /// Original enthalpy adjoint's total radiation-feedback scalar allowance.
    pub max_radiation_feedback_entries: usize,
}

/// Preserve the original typed refusal, including material corners and budgets.
#[derive(Debug)]
pub enum ScanSensitivityError {
    /// Primal/history/source assembly or phase admission refused.
    History(HeterogeneousEnthalpyError),
    /// Original checked implicit-function or transpose solver refused.
    Adjoint(EnthalpyAdjointError),
    /// Moving-source derivative, scalar admission, or cancellation refused.
    Source(ConductionError),
}
impl std::fmt::Display for ScanSensitivityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self { Self::History(e) => std::fmt::Display::fmt(e, f),
            Self::Adjoint(e) => std::fmt::Display::fmt(e, f), Self::Source(e) => std::fmt::Display::fmt(e, f) }
    }
}
impl std::error::Error for ScanSensitivityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::History(e) => Some(e), Self::Adjoint(e) => Some(e), Self::Source(e) => Some(e) }
    }
}
impl From<HeterogeneousEnthalpyError> for ScanSensitivityError {
    fn from(e: HeterogeneousEnthalpyError) -> Self { Self::History(e) }
}
impl From<EnthalpyAdjointError> for ScanSensitivityError {
    fn from(e: EnthalpyAdjointError) -> Self { Self::Adjoint(e) }
}
impl From<ConductionError> for ScanSensitivityError {
    fn from(e: ConductionError) -> Self { Self::Source(e) }
}

/// All sources share the SAME physical transpose; overlapping heaters add.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanStepGradient {
    /// Original source-density, history, reference-density and true-residual data.
    /// `transport.previous_specific_enthalpy` seeds the preceding fixed-grid step.
    pub storage: HeterogeneousEnthalpyStepGradient,
    /// Independent instantaneous controls, in source declaration order.
    pub sources: Vec<MovingHeatGradient>,
}

/// Primal report and immutable derivative binding for one accepted scan step.
/// Only the original mesh is borrowed; callers may retain a bounded collection
/// for a reverse sweep. The run itself never grows a trajectory tape.
#[derive(Debug)]
pub struct ScanStepLinearization<'m> {
    report: ScanStep,
    transport: HeterogeneousEnthalpyStepLinearization<'m>,
    sources: Vec<MovingHeatLinearization<'m>>,
    start_time_s: f64,
    end_time_s: f64,
}
impl ScanStepLinearization<'_> {
    /// Original full physical report, retaining source and radiation observations.
    #[must_use]
    pub const fn report(&self) -> &ScanStep { &self.report }
    /// Physical beginning of this fixed differentiation interval, seconds.
    #[must_use]
    pub const fn start_time_s(&self) -> f64 { self.start_time_s }
    /// Physical endpoint, held fixed by every returned derivative, seconds.
    #[must_use]
    pub const fn end_time_s(&self) -> f64 { self.end_time_s }
    /// Rechecked physical tangent, including any nonlinear radiation feedback.
    #[must_use]
    pub fn transport(&self) -> &HeterogeneousEnthalpyStepLinearization<'_> { &self.transport }

    /// Pull back an objective seed in nodal specific-enthalpy coordinates.
    /// Complete radiation feedback enters the common transpose, but reservoir
    /// and emissivity PARAMETER gradients are not returned by this source API.
    /// Geometry, phase charts, conductivity and boundary histories remain fixed.
    pub fn pullback(&self, cx: &Cx<'_>, seed: &[f64], linear: LinearConfig)
        -> Result<ScanStepGradient, ScanSensitivityError>
    {
        let storage = self.transport.pullback(cx, seed, linear)?;
        let mut sources = Vec::with_capacity(self.sources.len());
        for (i, binding) in self.sources.iter().enumerate() {
            poll(cx, i)?;
            sources.push(binding.pullback(cx, &storage.transport.source_density)?);
        }
        poll(cx, sources.len())?;
        Ok(ScanStepGradient { storage, sources })
    }

    /// Pull back a nodal temperature objective through the original dT/dh.
    /// Plateau interiors correctly have zero temperature sensitivity even when
    /// their enthalpy and liquid-fraction state can still change with heating.
    pub fn pullback_temperature(&self, cx: &Cx<'_>, seed: &[f64], linear: LinearConfig)
        -> Result<ScanStepGradient, ScanSensitivityError>
    {
        let enthalpy_seed = self.transport.transport().temperature_pullback(cx, seed)?;
        self.pullback(cx, &enthalpy_seed, linear)
    }
}

impl<'h, 'm, 'c> ScanningEnthalpyRun<'h, 'm, 'c> {
    /// Advance and bind derivatives at one explicit right-continuous endpoint.
    /// The original primal solver runs once. Publication occurs only after
    /// source binding and physical residual/energy/smoothness checks all pass.
    /// A valid primal at a phase kink may therefore refuse this derivative API;
    /// ordinary advance_to remains available and unchanged.
    ///
    /// # Errors
    /// Original history/motion/solver refusal, derivative memory admission,
    /// unsupported classical derivative, or cancellation. Every failure leaves
    /// the accepted caloric state, clock and step counter unchanged.
    pub fn advance_linearized_to(&mut self, cx: &Cx<'_>, time_s: f64,
        config: EnthalpyStepConfig, budget: ScanSensitivityBudget)
        -> Result<ScanStepLinearization<'m>, ScanSensitivityError>
    {
        self.admit(cx, time_s)?;
        if self.motion_endpoint(cx, self.run.time_s(), time_s)? < time_s {
            return Err(EnthalpyAdjointError::InvalidInput("linearized scan endpoint exceeds its selected motion limit").into());
        }
        let binding = self.linearize_endpoint(cx, time_s, false, config, budget)?;
        self.run.publish(cx, time_s, binding.report.solution())?;
        Ok(binding)
    }

    /// Advance to the next event/motion-limited endpoint and bind that interval.
    /// Held loads and their derivatives use the identical pre-switch values.
    /// Returned gradients hold the selected endpoint and event side FIXED; they
    /// are not derivatives of which event or motion limit selected the interval.
    pub fn advance_linearized_event_step(&mut self, cx: &Cx<'_>, target_s: f64,
        max_step_s: f64, config: EnthalpyStepConfig, budget: ScanSensitivityBudget)
        -> Result<ScanStepLinearization<'m>, ScanSensitivityError>
    {
        self.admit(cx, target_s)?;
        let start = self.run.time_s();
        let end = events::next_endpoint(cx, self.run.problem.boundary, &self.run.problem.source,
            self.signals(), start, target_s, max_step_s)?;
        let end = self.motion_endpoint(cx, start, end)?;
        let binding = self.linearize_endpoint(cx, end, true, config, budget)?;
        self.run.publish(cx, end, binding.report.solution())?;
        Ok(binding)
    }

    fn linearize_endpoint(&self, cx: &Cx<'_>, end: f64, left: bool,
        config: EnthalpyStepConfig, budget: ScanSensitivityBudget)
        -> Result<ScanStepLinearization<'m>, ScanSensitivityError>
    {
        poll(cx, 0)?;
        let mesh = self.run.problem.boundary.mesh();
        if self.sources.len().checked_mul(mesh.vertex_count())
            .is_none_or(|n| n > budget.max_source_vertices)
        {
            return Err(EnthalpyAdjointError::InvalidInput("scan source tangent storage exceeds its explicit source-vertex cap").into());
        }
        let mut sources = Vec::with_capacity(self.sources.len());
        for (i, source) in self.sources.iter().enumerate() {
            poll(cx, i)?;
            sources.push(if left { source.linearize_before(cx, mesh, end)? }
                else { source.linearize_at(cx, mesh, end)? });
        }
        let start = self.run.time_s();
        let report = self.solve(cx, &self.run.specific_h, start, end, left, config)?;
        // Immutable declarations and physical time are identical to the primal.
        // The existing binder independently reassembles the complete endpoint.
        let boundary = if left { events::boundary_before(cx, self.run.problem.boundary, end)? }
            else { self.run.problem.boundary.at(cx, end)? };
        natural_only(&boundary)?;
        let (source, _) = self.source_at(cx, end, left)?;
        let problem = ConductionProblem { mesh, boundary: &boundary, source: &source,
            material: self.run.problem.material, element_materials: self.run.problem.element_materials };
        let accepted = report.solution().clone();
        let transport = match &self.radiation {
            None => self.run.storage.linearize_accepted(cx, problem, self.run.problem.interfaces,
                &self.run.specific_h, end - start, config, accepted)?,
            Some(radiation) => {
                let patches = if left { patches_before(cx, &radiation.patches, end)? }
                    else {
                        let mut patches = Vec::with_capacity(radiation.patches.len());
                        for (i, patch) in radiation.patches.iter().enumerate() {
                            poll(cx, i)?;
                            patches.push(patch.at(end)?);
                        }
                        patches
                    };
                self.run.storage.linearize_accepted_storage_with_ambient_radiation(cx, problem,
                    self.run.problem.interfaces, &self.run.specific_h, end - start, config,
                    &patches, accepted, budget.max_radiation_feedback_entries)?
            }
        };
        poll(cx, self.sources.len())?;
        Ok(ScanStepLinearization { report, transport, sources, start_time_s: start, end_time_s: end })
    }
}
