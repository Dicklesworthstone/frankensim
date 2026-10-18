//! Advance scenario boundary/source histories through the existing solid stepper.
//!
//! All loads are sampled at the explicitly requested END time, including the
//! scenario owner's right-continuous held tables. This is backward-Euler
//! endpoint quadrature, not exact integration of a discontinuous schedule.
//! There is no hidden event splitting, timestep selection, continuous-peak
//! bound or Crank-Nicolson claim. The caller owns the time grid.
//!
//! The same old field feeds every attempted step. Prescribed-temperature
//! changes use BackwardEuler::advance_prescribed, including boundary capacity
//! in reaction heat. Failed sampling, cancellation, solve or energy checks
//! leave temperature, time and the accepted-step budget unchanged.

use super::{ThermalBoundaryHistory, poll};
use crate::{ConductionError, ConductionProblem, ConductivityModel, ElementMaterials,
    ScalarField, ThermalInterfaces, VOLUMETRIC_SOURCE_DIMS};
use crate::transient::VolumetricHeatCapacity;
use crate::transient::backward_euler::{BackwardEuler, StepConfig, StepSolution};
use fs_exec::Cx;
use fs_qty::Dims;
use fs_scenario::signal::TimeSignal;

/// Explicit volumetric generation, in W/m³, never an untyped total power.
#[derive(Debug, Clone, PartialEq)]
pub enum HeatSourceHistory {
    /// Original uniform or nodal field, unchanged over time.
    Fixed(ScalarField),
    /// Spatially uniform time history whose values have dimensions W/m³.
    Uniform(TimeSignal),
    /// A fixed W/m³ field multiplied by a dimensionless time history. This
    /// preserves a caller's nodal heater footprint while varying its workload.
    Scaled {
        /// Fixed spatial generation field, W/m³.
        field: ScalarField,
        /// Explicit dimensionless multiplier, including a declared zero/off state.
        scale: TimeSignal,
    },
}

impl HeatSourceHistory {
    fn at(&self, time_s: f64, vertices: usize) -> Result<ScalarField, ConductionError> {
        let field = match self {
            Self::Fixed(field) => field.clone(),
            Self::Uniform(signal) => ScalarField::Uniform(signal_value(signal, time_s, VOLUMETRIC_SOURCE_DIMS)?),
            Self::Scaled { field, scale } => {
                field.validate("heat source footprint", vertices)?;
                let factor = signal_value(scale, time_s, Dims([0; 6]))?;
                match field {
                    ScalarField::Uniform(value) => ScalarField::Uniform(value * factor),
                    ScalarField::Nodal(values) => ScalarField::Nodal(values.iter().map(|v| v * factor).collect()),
                }
            }
        };
        field.validate("history heat source", vertices)?;
        Ok(field)
    }
}

/// Fixed spatial/material data and the explicit histories that drive them.
/// No participant may be changed behind a running solve's immutable borrows.
pub struct ThermalHistoryProblem<'h, 'm> {
    /// Boundary histories, including the exact mesh and initial instant.
    pub boundary: &'h ThermalBoundaryHistory<'m>,
    /// Original conductivity model; must be temperature independent.
    pub material: &'h ConductivityModel,
    /// Optional exact per-element assignment; each model must be temperature independent.
    pub element_materials: Option<&'h ElementMaterials>,
    /// Optional matching-face thermal contacts, passed unchanged to the stepper.
    pub interfaces: Option<&'h ThermalInterfaces>,
    /// Fixed, uniform time-varying, or scaled-footprint generation.
    pub source: HeatSourceHistory,
}

/// In-memory, resumable thermal trajectory with an explicit lifetime step budget.
///
/// This is a driver of the existing linear backward-Euler owner, not a new
/// integrator. It retains only the accepted field; callers retain any desired
/// StepSolution history. Durable checkpointing and nonlinear/radiative transport
/// are not added by this interface.
pub struct ThermalHistoryRun<'h, 'm> {
    problem: ThermalHistoryProblem<'h, 'm>,
    stepper: BackwardEuler<'m>,
    temperature: Vec<f64>,
    time_s: f64,
    accepted_steps: usize,
    max_steps: usize,
}

impl<'h, 'm> ThermalHistoryRun<'h, 'm> {
    /// Start at the boundary history's explicit initial instant with uniform
    /// volumetric capacity. Initial temperatures must already match prescribed
    /// data: admission never resets boundary storage without accounting for it.
    ///
    /// # Errors
    /// Bad initial field/budget, incompatible initial prescribed temperatures,
    /// invalid source/history, nonlinear material, or capacity/cancellation refusal.
    pub fn new(
        cx: &Cx<'_>,
        problem: ThermalHistoryProblem<'h, 'm>,
        capacity: VolumetricHeatCapacity,
        initial: &[f64],
        max_steps: usize,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        if max_steps == 0 { return Err(invalid("the accepted-step budget must be positive")); }
        let mesh = problem.boundary.mesh();
        if initial.len() != mesh.vertex_count() {
            return Err(ConductionError::FieldLength {
                field: "initial history temperature", expected: mesh.vertex_count(), found: initial.len(),
            });
        }
        for (i, &value) in initial.iter().enumerate() {
            if i % 512 == 0 { poll(cx, i)?; }
            crate::require_finite("initial history temperature", value)?;
        }
        if let Some(materials) = problem.element_materials { materials.validate_for(mesh)?; }
        for e in 0..mesh.element_count() {
            if e % 512 == 0 { poll(cx, e)?; }
            let material = match problem.element_materials {
                Some(materials) => materials.model_for(e)?, None => problem.material,
            };
            if material.is_temperature_dependent() {
                return Err(invalid("history run requires temperature-independent conductivity; no frozen k(T) approximation is implicit"));
            }
        }
        let time_s = problem.boundary.initial_time_s();
        let boundary = problem.boundary.at(cx, time_s)?;
        for (i, &(vertex, prescribed)) in boundary.dirichlet().iter().enumerate() {
            if i % 512 == 0 { poll(cx, i)?; }
            if initial[vertex] != prescribed {
                return Err(invalid("initial field must match prescribed temperatures at the explicit initial time"));
            }
        }
        problem.source.at(time_s, mesh.vertex_count())?;
        let stepper = BackwardEuler::uniform(cx, mesh, capacity)?;
        poll(cx, 0)?;
        Ok(Self { problem, stepper, temperature: initial.to_vec(), time_s, accepted_steps: 0, max_steps })
    }

    /// Accepted full nodal field, including prescribed vertices, in kelvin.
    #[must_use]
    pub fn temperature(&self) -> &[f64] { &self.temperature }

    /// Time of the accepted field, in seconds.
    #[must_use]
    pub const fn time_s(&self) -> f64 { self.time_s }

    /// Number of successfully published steps; refused attempts do not count.
    #[must_use]
    pub const fn accepted_steps(&self) -> usize { self.accepted_steps }

    /// Explicit total accepted-step allowance, not a per-call budget.
    #[must_use]
    pub const fn max_steps(&self) -> usize { self.max_steps }

    /// Increase the lifetime budget without resetting temperature or physical time.
    ///
    /// # Errors
    /// A total allowance that does not exceed the current one.
    pub fn extend_step_budget(&mut self, max_steps: usize) -> Result<(), ConductionError> {
        if max_steps <= self.max_steps { return Err(invalid("step-budget extension must increase the total allowance")); }
        self.max_steps = max_steps;
        Ok(())
    }

    /// Advance to one explicitly requested, strictly later finite endpoint.
    ///
    /// The returned production StepSolution includes source/Neumann/Robin heat,
    /// prescribed-temperature reaction, storage change and independent energy
    /// closure. No separate audit formula or widened solve tolerance is used.
    /// Resume after cancellation with a fresh Cx and the same endpoint; after
    /// exhaustion, explicitly extend the budget. Result errors publish nothing.
    ///
    /// # Errors
    /// Exhaustion, invalid time/data, cancellation, or any production stepper
    /// refusal, including failure of the declared residual/energy gates.
    pub fn advance_to(
        &mut self,
        cx: &Cx<'_>,
        time_s: f64,
        config: StepConfig,
    ) -> Result<StepSolution, ConductionError> {
        poll(cx, self.accepted_steps)?;
        if self.accepted_steps >= self.max_steps { return Err(invalid("thermal history step budget exhausted")); }
        let dt_s = time_s - self.time_s;
        if !time_s.is_finite() || !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(invalid("thermal history endpoint must be finite and strictly later with a representable positive interval"));
        }
        let boundary = self.problem.boundary.at(cx, time_s)?;
        let mesh = self.problem.boundary.mesh();
        let source = self.problem.source.at(time_s, mesh.vertex_count())?;
        poll(cx, self.accepted_steps)?;
        let solution = self.stepper.advance_prescribed(cx, ConductionProblem {
            mesh, boundary: &boundary, material: self.problem.material, source: &source,
            element_materials: self.problem.element_materials,
        }, self.problem.interfaces, &self.temperature, dt_s, config)?;
        // Prepare the published copy before changing any accepted state.
        let temperature = solution.temperature.clone();
        poll(cx, self.accepted_steps)?;
        self.temperature = temperature;
        self.time_s = time_s;
        self.accepted_steps += 1;
        Ok(solution)
    }
}

fn signal_value(signal: &TimeSignal, time_s: f64, dims: Dims) -> Result<f64, ConductionError> {
    let value = signal.eval(time_s).map_err(|error| invalid(&format!("heat-source history: {error}")))?;
    if value.dims != dims {
        return Err(ConductionError::Dimensions {
            context: "heat-source history".into(), expected: dims.0, found: value.dims.0,
        });
    }
    Ok(value.value)
}

fn invalid(what: &str) -> ConductionError {
    ConductionError::Config { parameter: "thermal history run", what: what.into() }
}
