//! Advance scenario boundary/source histories through the existing solid steppers.
//!
//! All loads are sampled at the explicitly requested END time, including the
//! scenario owner's right-continuous held tables. This is backward-Euler
//! endpoint quadrature, not exact integration of a discontinuous schedule.
//! The original advance_to methods keep that endpoint rule. The opt-in
//! advance_event_step methods resolve schedule knots and use their left limits.
//! Neither path supplies a temporal error bound or a continuous-peak bound.
//!
//! The same old field feeds every attempted step. Prescribed-temperature
//! changes use the existing endpoint solvers, including boundary capacity
//! in reaction heat. Failed sampling, cancellation, solve or energy checks
//! leave temperature, time and the accepted-step budget unchanged.

/// Sourced radiative surfaces and independent reservoir-temperature histories.
pub mod radiation;

/// Scheduled reference-mass enthalpy transport, including latent phase changes.
pub mod enthalpy;

mod capacity;
pub use capacity::ThermalHistoryCapacity;

/// Event-aligned steps that do not skip scheduled loads or backdate their jumps.
pub mod events;

/// Moving power-normalized source fields for the existing thermal steppers.
pub mod moving;

use super::{ThermalBoundaryHistory, poll};
use crate::{ConductionError, ConductionProblem, ConductivityModel, ElementMaterials,
    ScalarField, ThermalBoundary, ThermalInterfaces, VOLUMETRIC_SOURCE_DIMS};
use crate::transient::VolumetricHeatCapacity;
use crate::transient::backward_euler::{BackwardEuler, NonlinearStepConfig,
    NonlinearStepSolution, StepConfig, StepSolution};
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
    /// Original conductivity model. k(T) requires an explicit nonlinear method.
    pub material: &'h ConductivityModel,
    /// Optional exact per-element assignment, subject to the same method policy.
    pub element_materials: Option<&'h ElementMaterials>,
    /// Optional matching-face thermal contacts, passed unchanged to the stepper.
    pub interfaces: Option<&'h ThermalInterfaces>,
    /// Fixed, uniform time-varying, or scaled-footprint generation.
    pub source: HeatSourceHistory,
}

/// Immutable choice of endpoint transport and its nonlinear work budget.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalHistoryMethod {
    /// Constant conductivity only; retains the original linear solve and refusal.
    Linear,
    /// Newton/FGMRES on the actual endpoint k(T), not a frozen old-temperature law.
    Nonlinear(NonlinearStepConfig),
}

/// Original solver response, including nonlinear evidence when that path is used.
#[derive(Debug, Clone)]
pub enum ThermalHistoryStep {
    /// Constant-conductivity endpoint.
    Linear(StepSolution),
    /// Nonlinear endpoint with its actual Newton work and residual threshold.
    Nonlinear(NonlinearStepSolution),
}

impl ThermalHistoryStep {
    /// Accepted temperatures, fluxes, storage change and discrete energy closure.
    #[must_use]
    pub fn solution(&self) -> &StepSolution {
        match self { Self::Linear(step) => step, Self::Nonlinear(solved) => &solved.step }
    }

    /// Consume the response when only the common physical observations are needed.
    #[must_use]
    pub fn into_solution(self) -> StepSolution {
        match self { Self::Linear(step) => step, Self::Nonlinear(solved) => solved.step }
    }
}

/// In-memory, resumable thermal trajectory with an explicit lifetime step budget.
///
/// This drives the existing backward-Euler owners, not a new integrator. It
/// retains only the accepted field; callers retain any desired step reports.
/// Heat capacity and contact resistance remain temperature independent.
/// No durable checkpointing or latent-heat model is added by this interface.
pub struct ThermalHistoryRun<'h, 'm> {
    problem: ThermalHistoryProblem<'h, 'm>,
    stepper: BackwardEuler<'m>,
    method: ThermalHistoryMethod,
    temperature: Vec<f64>,
    time_s: f64,
    accepted_steps: usize,
    max_steps: usize,
}

impl<'h, 'm> ThermalHistoryRun<'h, 'm> {
    /// Start a LINEAR run at the boundary history's explicit initial instant.
    /// Initial temperatures must already match prescribed data: admission never
    /// resets boundary storage without accounting for it. k(T) remains refused;
    /// use [`Self::new_with_method`] to select a nonlinear policy explicitly.
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
        Self::new_with_method(cx, problem, capacity, initial, max_steps, ThermalHistoryMethod::Linear)
    }

    /// Start with an explicitly selected linear or nonlinear endpoint method.
    ///
    /// Nonlinear transport evaluates the original conductivity and its tangent
    /// at every candidate field, including NEW prescribed temperatures. The
    /// Newton policy remains fixed for this run; each step's `StepConfig` retains
    /// the production total Krylov budget and independent joule-closure gate.
    /// Material range checks stay with the solver; no curve is extrapolated here.
    ///
    /// # Errors
    /// Invalid nonlinear policy, or any admission refusal from [`Self::new`]
    /// except that the nonlinear method admits temperature-dependent conductivity.
    pub fn new_with_method(
        cx: &Cx<'_>,
        problem: ThermalHistoryProblem<'h, 'm>,
        capacity: VolumetricHeatCapacity,
        initial: &[f64],
        max_steps: usize,
        method: ThermalHistoryMethod,
    ) -> Result<Self, ConductionError> {
        Self::new_with_capacity(cx, problem, ThermalHistoryCapacity::Uniform(capacity),
            initial, max_steps, method)
    }

    /// Start with uniform, element-wise, or material-ID-bound heat storage.
    ///
    /// Every capacity is a declared or source-resolved constant rho*c_p in J/(m^3 K). The original
    /// lumped P1 assembly is reused; values are never averaged across materials
    /// before assembly. Material-ID binding uses the SAME element assignment
    /// as conductivity and refuses missing, duplicate or unused capacity rows.
    /// The resulting capacity matrix stays fixed through nonlinear iterations,
    /// radiation, cancellation and resume. No temperature-dependent capacity or
    /// latent-heat approximation is introduced.
    ///
    /// # Errors
    /// The same refusals as [`Self::new_with_method`], plus mismatched element
    /// counts, invalid material bindings or unrepresentable capacity assembly.
    pub fn new_with_capacity(
        cx: &Cx<'_>,
        problem: ThermalHistoryProblem<'h, 'm>,
        capacity: ThermalHistoryCapacity<'_>,
        initial: &[f64],
        max_steps: usize,
        method: ThermalHistoryMethod,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        if let ThermalHistoryMethod::Nonlinear(policy) = method { policy.validate()?; }
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
            if matches!(method, ThermalHistoryMethod::Linear) && material.is_temperature_dependent() {
                return Err(invalid("history run requires an explicit nonlinear method for k(T); no frozen-conductivity approximation is implicit"));
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
        let stepper = capacity.prepare(cx, mesh, problem.element_materials)?;
        stepper.validate_capacity_temperature(cx, initial)?;
        poll(cx, 0)?;
        Ok(Self { problem, stepper, method, temperature: initial.to_vec(), time_s, accepted_steps: 0, max_steps })
    }

    /// Accepted full nodal field, including prescribed vertices, in kelvin.
    #[must_use]
    pub fn temperature(&self) -> &[f64] { &self.temperature }

    /// Original card identities, density/Cp claims and admitted temperature spans.
    /// Empty for caller-declared capacities.
    #[must_use]
    pub fn material_capacities(&self) -> &[crate::transient::MaterialHeatCapacity] {
        self.stepper.material_capacities()
    }

    /// Time of the accepted field, in seconds.
    #[must_use]
    pub const fn time_s(&self) -> f64 { self.time_s }

    /// Number of successfully published steps; refused attempts do not count.
    #[must_use]
    pub const fn accepted_steps(&self) -> usize { self.accepted_steps }

    /// Explicit total accepted-step allowance, not a per-call budget.
    #[must_use]
    pub const fn max_steps(&self) -> usize { self.max_steps }

    /// Original immutable endpoint method and nonlinear policy.
    #[must_use]
    pub const fn method(&self) -> ThermalHistoryMethod { self.method }

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
    /// The common response retains physical flux/storage/energy observations.
    /// Use [`Self::advance_to_report`] to also retain the nonlinear solver's
    /// Newton/backtrack counts and its independently checked residual gate.
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
        self.advance_to_report(cx, time_s, config).map(ThermalHistoryStep::into_solution)
    }

    /// Advance with the full response from the explicitly selected solver.
    ///
    /// Nonlinear iterations all reuse the SAME old field and endpoint loads.
    /// Only a residual- and energy-admitted endpoint advances physical time.
    /// Constant materials may use the nonlinear method for an independent
    /// comparison, but are not silently redirected to the linear implementation.
    ///
    /// # Errors
    /// The same transactional refusals as [`Self::advance_to`].
    pub fn advance_to_report(
        &mut self, cx: &Cx<'_>, time_s: f64, config: StepConfig,
    ) -> Result<ThermalHistoryStep, ConductionError> {
        let (dt_s, boundary, source) = self.prepare_endpoint(cx, time_s)?;
        let problem = ConductionProblem {
            mesh: self.problem.boundary.mesh(), boundary: &boundary,
            material: self.problem.material, source: &source,
            element_materials: self.problem.element_materials,
        };
        let report = match self.method {
            ThermalHistoryMethod::Linear => ThermalHistoryStep::Linear(
                self.stepper.advance_prescribed(cx, problem, self.problem.interfaces,
                    &self.temperature, dt_s, config)?),
            ThermalHistoryMethod::Nonlinear(policy) => ThermalHistoryStep::Nonlinear(
                self.stepper.advance_nonlinear_prescribed(cx, problem, self.problem.interfaces,
                    &self.temperature, dt_s, config, policy)?),
        };
        self.publish(cx, time_s, report.solution())?;
        Ok(report)
    }

    fn prepare_endpoint(&self, cx: &Cx<'_>, time_s: f64)
        -> Result<(f64, ThermalBoundary, ScalarField), ConductionError>
    {
        poll(cx, self.accepted_steps)?;
        if self.accepted_steps >= self.max_steps { return Err(invalid("thermal history step budget exhausted")); }
        let dt_s = time_s - self.time_s;
        if !time_s.is_finite() || !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(invalid("thermal history endpoint must be finite and strictly later with a representable positive interval"));
        }
        let boundary = self.problem.boundary.at(cx, time_s)?;
        let source = self.problem.source.at(time_s, self.problem.boundary.mesh().vertex_count())?;
        poll(cx, self.accepted_steps)?;
        Ok((dt_s, boundary, source))
    }

    fn publish(&mut self, cx: &Cx<'_>, time_s: f64, solution: &StepSolution) -> Result<(), ConductionError> {
        // Prepare the published copy before changing any accepted state.
        let temperature = solution.temperature.clone();
        poll(cx, self.accepted_steps)?;
        self.temperature = temperature;
        self.time_s = time_s;
        self.accepted_steps += 1;
        Ok(())
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
