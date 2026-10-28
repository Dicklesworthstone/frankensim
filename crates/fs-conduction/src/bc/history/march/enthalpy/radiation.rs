//! Radiation-driven phase change with immutable material and reservoir histories.
//!
//! The existing heterogeneous enthalpy/radiation solver owns the full implicit
//! feedback, reference-mass storage, phase-dependent conductivity and heat split.
//! The wrapper samples original reservoir histories and publishes only its
//! energy-admitted enthalpy/temperature/fraction endpoint. No inner iteration
//! advances physical time, and no latent heat is deposited a second time.

mod adaptive;

use super::{EnthalpyHistoryRun, conduction, invalid, poll};
use super::super::radiation::{RadiationRegionHistory, admit_histories, patches_before};
use super::super::events;
use crate::{AmbientRadiationConfig, ConductionProblem};
use crate::transient::enthalpy::EnthalpyStepConfig;
use crate::transient::enthalpy::heterogeneous::HeterogeneousEnthalpyError;
use crate::transient::enthalpy::radiation::EnthalpyRadiationStepSolution;
use fs_exec::Cx;

/// Persistent radiative phase-change trajectory. No mutable escape can advance
/// the underlying enthalpy without its attached surface/reservoir declarations.
/// The model remains the original area-mean gray patch facing a black reservoir,
/// not an enclosure or a pointwise fourth-power surface quadrature.
pub struct RadiativeEnthalpyHistoryRun<'h, 'm, 'c> {
    run: EnthalpyHistoryRun<'h, 'm, 'c>,
    patches: Vec<RadiationRegionHistory>,
    config: AmbientRadiationConfig,
}

impl<'h, 'm, 'c> RadiativeEnthalpyHistoryRun<'h, 'm, 'c> {
    /// Attach sourced emissivities and independent Kelvin-valued reservoirs
    /// before the first accepted step. Both Robin and explicit zero-flux traces
    /// use the temperature runner's SAME binding rules. Radiation-only traces
    /// add no fictitious convection; material validity remains solver-owned.
    ///
    /// # Errors
    /// Invalid patch/time/outer policy, duplicate or unknown regions, an advanced
    /// run, or cancellation. No physical state is advanced by construction.
    pub fn new(cx: &Cx<'_>, run: EnthalpyHistoryRun<'h, 'm, 'c>,
        patches: Vec<RadiationRegionHistory>, config: AmbientRadiationConfig)
        -> Result<Self, HeterogeneousEnthalpyError>
    {
        poll(cx, 0).map_err(conduction)?;
        config.validate().map_err(conduction)?;
        if run.accepted_steps() != 0 {
            return Err(invalid("radiation must be attached before the first enthalpy step"));
        }
        admit_histories(cx, run.problem.boundary, run.time_s(), &patches).map_err(conduction)?;
        Ok(Self { run, patches, config })
    }

    /// Accepted enthalpy, temperature, liquid fractions, clock and material state.
    #[must_use]
    pub fn conduction(&self) -> &EnthalpyHistoryRun<'h, 'm, 'c> { &self.run }

    /// Original emissivity receipts and independent reservoir schedules.
    #[must_use]
    pub fn patches(&self) -> &[RadiationRegionHistory] { &self.patches }

    /// Unchanged radiative iteration policy, also retained across retries.
    #[must_use]
    pub const fn config(&self) -> AmbientRadiationConfig { self.config }

    /// Extend the physical-step allowance without dropping radiation or phase.
    ///
    /// # Errors
    /// The same non-increasing-budget refusal as the underlying run.
    pub fn extend_step_budget(&mut self, max_steps: usize) -> Result<(), HeterogeneousEnthalpyError> {
        self.run.extend_step_budget(max_steps)
    }

    /// Advance all loads and reservoirs at one explicitly requested endpoint.
    ///
    /// Returns the original solver's complete response, including separate
    /// convection, nonlinear radiation, phase fields, cumulative Newton/Krylov
    /// work and the independently reassembled physical residual/energy checks.
    /// The frozen inner Robin report is not relabeled as pure convection.
    ///
    /// # Errors
    /// Original endpoint, budget, sampling, phase, solver or cancellation refusal.
    /// None changes accepted enthalpy, temperature, liquid fraction or time.
    pub fn advance_to(&mut self, cx: &Cx<'_>, time_s: f64, config: EnthalpyStepConfig)
        -> Result<EnthalpyRadiationStepSolution, HeterogeneousEnthalpyError>
    {
        let (dt_s, boundary, source) = self.run.prepare_endpoint(cx, time_s)?;
        let mut patches = Vec::with_capacity(self.patches.len());
        for (at, history) in self.patches.iter().enumerate() {
            poll(cx, at).map_err(conduction)?;
            patches.push(history.at(time_s).map_err(conduction)?);
        }
        let solved = self.run.storage.advance_with_ambient_radiation(cx, ConductionProblem {
            mesh: self.run.problem.boundary.mesh(), boundary: &boundary,
            material: self.run.problem.material, element_materials: self.run.problem.element_materials,
            source: &source,
        }, self.run.problem.interfaces, &self.run.specific_h, dt_s, config, &patches, self.config)?;
        self.run.publish(cx, time_s, &solved.conduction)?;
        Ok(solved)
    }

    /// Advance ONE event-aligned phase/radiation step, preserving every source.
    /// Reservoir, heat-source and natural-boundary knots all constrain the next
    /// endpoint. Held values use the exact left limit of this interval; latent
    /// heat and phase remain in the original constitutive charts and solver.
    /// Repeat until `conduction().time_s() == target_s` to traverse a window.
    ///
    /// # Errors
    /// Invalid target/cap or the existing budget, domain, solve and cancellation
    /// refusals. A failed interval publishes none of its caloric fields or time.
    pub fn advance_event_step(&mut self, cx: &Cx<'_>, target_s: f64,
        max_step_s: f64, config: EnthalpyStepConfig)
        -> Result<events::HistoryEventStep<EnthalpyRadiationStepSolution>, HeterogeneousEnthalpyError>
    {
        poll(cx, self.run.accepted_steps()).map_err(conduction)?;
        if self.run.accepted_steps() >= self.run.max_steps() {
            return Err(invalid("enthalpy history step budget exhausted"));
        }
        let start_time_s = self.run.time_s();
        let end_time_s = events::next_endpoint(cx, self.run.problem.boundary, &self.run.problem.source,
            self.patches.iter().map(|p| &p.ambient_temperature), start_time_s, target_s, max_step_s)
            .map_err(conduction)?;
        let boundary = events::boundary_before(cx, self.run.problem.boundary, end_time_s).map_err(conduction)?;
        self.run.storage.admit_prescribed_boundary(cx, &boundary)?;
        let source = events::source_before(&self.run.problem.source, end_time_s, self.run.specific_h.len())
            .map_err(conduction)?;
        let patches = patches_before(cx, &self.patches, end_time_s).map_err(conduction)?;
        let report = self.run.storage.advance_with_ambient_radiation(cx, ConductionProblem {
            mesh: self.run.problem.boundary.mesh(), boundary: &boundary,
            material: self.run.problem.material, element_materials: self.run.problem.element_materials,
            source: &source,
        }, self.run.problem.interfaces, &self.run.specific_h, end_time_s-start_time_s,
            config, &patches, self.config)?;
        self.run.publish(cx, end_time_s, &report.conduction)?;
        Ok(events::HistoryEventStep { start_time_s, end_time_s, report })
    }
}
