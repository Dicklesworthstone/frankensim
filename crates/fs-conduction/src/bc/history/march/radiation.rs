//! Scheduled ambient radiation coupled to the history runner's actual endpoint.
//!
//! Each patch keeps its material-card emissivity and an independent Kelvin-valued
//! black-reservoir history. The existing area-mean T^4 model and nonlinear
//! endpoint solver own the physics, material validity, heat split and energy
//! checks. Reservoir temperatures are not inferred from convection references.
//!
//! The wrapper owns an unadvanced run and immutable radiation declarations, so
//! later steps and resumes cannot accidentally omit or retune a radiating patch.
//! All schedules use explicit endpoint quadrature. This adds no enclosure,
//! occlusion, pointwise T(x)^4 integral, evolving emissivity or latent-heat model.

use super::{ThermalBoundaryHistory, ThermalHistoryMethod, ThermalHistoryRun, invalid, poll};
use crate::{AmbientRadiationConfig, AmbientRadiationPatch, ConductionError,
    ConductionProblem, ScalarField, SurfaceEmissivity, ThermalBc, TEMPERATURE_DIMS};
use crate::transient::backward_euler::{RadiationStepSolution, StepConfig};
use fs_exec::Cx;
use fs_scenario::signal::TimeSignal;
use std::collections::BTreeSet;

/// One sourced surface and its explicitly scheduled radiative surroundings.
#[derive(Debug, Clone, PartialEq)]
pub struct RadiationRegionHistory {
    /// Exact existing natural-boundary region, distinct for each radiating patch.
    pub region: String,
    /// Original material-card value, selection receipt and surface validity range.
    pub emissivity: SurfaceEmissivity,
    /// Black-reservoir absolute temperature, K; independent of convection air.
    pub ambient_temperature: TimeSignal,
}

impl RadiationRegionHistory {
    /// Sample the reservoir through the original signal and patch constructors.
    ///
    /// # Errors
    /// Malformed/nonfinite signal, wrong quantity dimensions, negative absolute
    /// temperature, invalid patch identity or an unrepresentable radiative law.
    /// Surface-material validity is additionally checked by the endpoint solver.
    pub fn at(&self, time_s: f64) -> Result<AmbientRadiationPatch, ConductionError> {
        let value = self.ambient_temperature.eval(time_s).map_err(|error|
            refused(&self.region, &format!("radiative reservoir history: {error}")))?;
        if value.dims != TEMPERATURE_DIMS {
            return Err(ConductionError::Dimensions {
                context: format!("radiative reservoir on {:?}", self.region),
                expected: TEMPERATURE_DIMS.0, found: value.dims.0,
            });
        }
        AmbientRadiationPatch::new(self.region.clone(), self.emissivity.clone(), value.value)
    }
}

/// A history run whose radiation model remains attached for its whole trajectory.
/// Conduction may be linear or explicitly nonlinear, as selected on the inner run.
pub struct RadiativeThermalHistoryRun<'h, 'm> {
    run: ThermalHistoryRun<'h, 'm>,
    patches: Vec<RadiationRegionHistory>,
    config: AmbientRadiationConfig,
}

impl<'h, 'm> RadiativeThermalHistoryRun<'h, 'm> {
    /// Attach radiation before the first accepted physical step.
    ///
    /// Each surface must own a nonempty, unique, uniform Robin or explicitly
    /// zero-flux trace. A zero-flux trace means no convection: no artificial h
    /// or convective reference is introduced. Initial declarations are checked
    /// here; changing histories are checked again at every actual endpoint.
    ///
    /// # Errors
    /// Already advanced run, empty/duplicate/unknown/unsupported patches,
    /// invalid histories or outer policy, and cancellation. No step is taken.
    pub fn new(
        cx: &Cx<'_>, run: ThermalHistoryRun<'h, 'm>,
        patches: Vec<RadiationRegionHistory>, config: AmbientRadiationConfig,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        config.validate()?;
        if run.accepted_steps() != 0 {
            return Err(invalid("radiation histories must be attached before the first accepted step"));
        }
        admit_histories(cx, run.problem.boundary, run.time_s(), &patches)?;
        Ok(Self { run, patches, config })
    }

    /// Read-only conduction state and immutable method; no mutable escape can
    /// advance the inner run without its attached radiation declarations.
    #[must_use]
    pub fn conduction(&self) -> &ThermalHistoryRun<'h, 'm> { &self.run }

    /// Original sourced surface and reservoir histories, unchanged by stepping.
    #[must_use]
    pub fn patches(&self) -> &[RadiationRegionHistory] { &self.patches }

    /// Explicit outer iteration policy, unchanged by retries and resume.
    #[must_use]
    pub const fn config(&self) -> AmbientRadiationConfig { self.config }

    /// Increase the total physical step budget without dropping radiation.
    ///
    /// # Errors
    /// A non-increase, with the same semantics as the inner run.
    pub fn extend_step_budget(&mut self, max_steps: usize) -> Result<(), ConductionError> {
        self.run.extend_step_budget(max_steps)
    }

    /// Advance all source, boundary and reservoir histories to one endpoint.
    ///
    /// Every radiative/Newton trial starts from the SAME accepted old field.
    /// The complete production response separates original convection from
    /// radiation and retains physical residual, reaction and joule-closure
    /// checks, material receipts and cumulative inner work. No balance budget
    /// is widened and no frozen inner report is relabeled as a physical one.
    ///
    /// # Errors
    /// History/domain/endpoint/budget/cancellation or existing radiative solver
    /// refusal. None changes accepted temperatures, time or the step counter.
    pub fn advance_to(&mut self, cx: &Cx<'_>, time_s: f64, config: StepConfig)
        -> Result<RadiationStepSolution, ConductionError>
    {
        let (dt_s, boundary, source) = self.run.prepare_endpoint(cx, time_s)?;
        let mut patches = Vec::with_capacity(self.patches.len());
        for (at, history) in self.patches.iter().enumerate() {
            poll(cx, at)?;
            patches.push(history.at(time_s)?);
        }
        let nonlinear = match self.run.method {
            ThermalHistoryMethod::Linear => None,
            ThermalHistoryMethod::Nonlinear(policy) => Some(policy),
        };
        let solution = self.run.stepper.advance_prescribed_with_ambient_radiation(cx,
            ConductionProblem { mesh: self.run.problem.boundary.mesh(), boundary: &boundary,
                material: self.run.problem.material, source: &source,
                element_materials: self.run.problem.element_materials },
            self.run.problem.interfaces, &self.run.temperature, dt_s, config,
            nonlinear, &patches, self.config)?;
        self.run.publish(cx, time_s, &solution.conduction)?;
        Ok(solution)
    }
}

// One binding contract for both temperature-capacity and primary-enthalpy runs.
pub(super) fn admit_histories(
    cx: &Cx<'_>, history: &ThermalBoundaryHistory<'_>, time_s: f64,
    patches: &[RadiationRegionHistory],
) -> Result<(), ConductionError> {
    poll(cx, 0)?;
    if patches.is_empty() {
        return Err(invalid("radiative history requires at least one declared patch; use the original run for no radiation"));
    }
    let boundary = history.at(cx, time_s)?;
    let mut seen = BTreeSet::new();
    for (at, patch) in patches.iter().enumerate() {
        poll(cx, at)?;
        let Some(index) = boundary.region_names().iter().position(|name| name == &patch.region) else {
            return Err(refused(&patch.region, "radiative patch has no matching boundary region"));
        };
        if !seen.insert(index) {
            return Err(refused(&patch.region, "radiative region is declared more than once"));
        }
        let admitted = match &boundary.conditions()[index] {
            ThermalBc::Robin { htc: ScalarField::Uniform(_), t_ref: ScalarField::Uniform(_) } => true,
            ThermalBc::Neumann { outward_flux: ScalarField::Uniform(flux) } => *flux == 0.0,
            _ => false,
        };
        if !admitted {
            return Err(refused(&patch.region, "radiation requires a uniform Robin or explicit zero-flux region"));
        }
        let mut owns_face = false;
        for slot in 0..history.mesh().boundary().len() {
            if slot % 512 == 0 { poll(cx, slot)?; }
            owns_face |= boundary.region_for(slot) == Some(index);
        }
        if !owns_face {
            return Err(refused(&patch.region, "radiative region owns no boundary faces"));
        }
        patch.at(time_s)?;
    }
    poll(cx, patches.len())
}

fn refused(region: &str, what: &str) -> ConductionError {
    ConductionError::Radiation { surface: region.into(), what: what.into(),
        fix: "bind a unique natural boundary to a sourced emissivity and a valid Kelvin reservoir history".into() }
}
