//! Spatial temperature reconstruction using the production P1 conduction solve.
//!
//! The base policy takes one backward-Euler step per interval; `substeps`
//! separates numerical resolution from observation times. Both use consistent heat
//! capacity, boundary loads, optional matching contacts and optional k(T).
//! Only free nodal temperatures are controls: prescribed values are lifted in
//! the primal and have zero tangent/cotangent, never optimized or penalized as
//! model error. Use the existing variational window and L-BFGS studies.
//!
//! Every forecast passes the conduction residual/energy gates. Inferred knot
//! corrections are statistical model discrepancies, not certified heat inputs;
//! the corrected history itself need not conserve energy. Geometry, capacity
//! laws, boundary partition and prescribed values are fixed across the window.
//! No air-network feedback, moving boundary or radiation derivative is implied.

use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, DofMap,
    ThermalBoundary, ThermalInterfaces};
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig,
    StepConfig, StepLinearization};
use fs_exec::Cx;
use fs_time::adaptive::adjoint::trajectory::TrajectoryGradient;
use crate::transient::variational::{WindowError, intervals::{IntervalScheme, IntervalTape}};

/// Bounded multistep forecasts with checkpointed production adjoints.
pub mod substeps;

/// Exact endpoint context shared by a substep's forward and derivative calls.
/// `interval` indexes observation knots; `index` indexes numerical steps inside
/// that interval. Data is evaluated at the endpoint, as required by backward
/// Euler. A discontinuous schedule must align with this explicitly fixed grid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConductionSubstep {
    pub interval: usize,
    pub index: usize,
    pub start: f64,
    pub end: f64,
}

/// Immutable model data for a parameter point. The capacity object and all
/// problems must use the policy's EXACT mesh. `problem(k)` supplies endpoint
/// material/source/boundary data for interval k; no time interpolation occurs.
/// The default parameter action declares these physical inputs independent of
/// parameters (which may still enter a WindowObjective as sensor parameters).
///
/// A parameterized model overrides parameter_pullback, for example contracting
/// step.source_density_pullback with a nodal source derivative, or using
/// step.capacity_multiplier_pullback for a log-capacity coordinate. Derivatives
/// must use the SAME coordinates as the parameter factory. Model and derivative
/// callbacks must remain unchanged throughout a tape's lifetime and retries.
pub trait ConductionWindowModel {
    fn engine(&self) -> &BackwardEuler<'_>;
    fn problem(&self, interval: usize) -> Result<ConductionProblem<'_>, ConductionError>;
    /// Endpoint data for a refined numerical step. The default freezes the
    /// original interval data. Override for a prescribed source/Robin schedule;
    /// a changing Dirichlet lift remains unsupported. Replay uses identical times.
    fn problem_at(&self, step: ConductionSubstep) -> Result<ConductionProblem<'_>, ConductionError> {
        self.problem(step.interval)
    }
    fn interfaces(&self) -> Option<&ThermalInterfaces> { None }
    fn parameter_count(&self) -> usize { 0 }
    fn parameter_pullback(&self, _interval: usize, _cx: &Cx<'_>,
        _step: &StepLinearization<'_>, _nodal_load_adjoint: &[f64],
        parameters: &mut [f64]) -> Result<(), ConductionError>
    {
        parameters.fill(0.0); Ok(())
    }
    /// Parameter action for exactly the endpoint data supplied by `problem_at`.
    /// Override together with that method when a schedule is parameterized.
    fn parameter_pullback_at(&self, time: ConductionSubstep, cx: &Cx<'_>,
        step: &StepLinearization<'_>, nodal_load_adjoint: &[f64],
        parameters: &mut [f64]) -> Result<(), ConductionError>
    {
        self.parameter_pullback(time.interval, cx, step, nodal_load_adjoint, parameters)
    }
}

/// Mesh/clock caps bound the amount of admitted assembly work. Production
/// conduction owns its sparse matrices, factor/preconditioner and Krylov memory;
/// the variational WindowControl separately bounds its control/result vectors.
/// These are not byte-exact whole-process allocation or latency guarantees.
#[derive(Debug, Clone)]
pub struct ConductionWindowConfig {
    pub step: StepConfig,
    pub nonlinear: Option<NonlinearStepConfig>,
    pub max_vertices: usize,
    pub max_elements: usize,
    pub max_intervals: usize,
    pub max_parameters: usize,
}

/// One fixed mesh, boundary lift, timetable and physical execution context.
/// The supplied Cx is used INSIDE all assembly/solve operations. The additional
/// window cancellation callback is polled between these operations; callers
/// needing cancellation within a solve must also request it through Cx's gate.
#[derive(Clone)]
pub struct ConductionWindowPolicy<'a, 'cx> {
    cx: &'a Cx<'cx>,
    mesh: &'a ConductionMesh,
    dofs: DofMap,
    times: Vec<f64>,
    config: ConductionWindowConfig,
}
fn poll(check: &mut dyn FnMut() -> bool) -> Result<(), WindowError> {
    if check() { Err(WindowError::Cancelled) } else { Ok(()) }
}
fn refusal(interval: usize, phase: &'static str, error: ConductionError) -> WindowError {
    match error {
        ConductionError::Cancelled { .. } => WindowError::Cancelled,
        other => WindowError::Integrator { interval, phase, diagnostic: other.to_string() },
    }
}
fn finite(values: &[f64]) -> bool { values.iter().all(|x| x.is_finite()) }
fn zeros(n: usize) -> Result<Vec<f64>, WindowError> {
    let mut values = Vec::new();
    values.try_reserve_exact(n).map_err(|_| WindowError::Allocation)?;
    values.resize(n, 0.0); Ok(values)
}
impl<'a, 'cx> ConductionWindowPolicy<'a, 'cx> {
    pub fn new(cx: &'a Cx<'cx>, mesh: &'a ConductionMesh, boundary: &ThermalBoundary,
        times: &[f64], config: ConductionWindowConfig, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, WindowError> {
        poll(cancelled)?;
        let n = mesh.vertex_count();
        if n == 0 || n > config.max_vertices || mesh.element_count() > config.max_elements
            || times.len() < 2 || times.len()-1 > config.max_intervals
            || !finite(times) || times.windows(2).any(|t| t[1] <= t[0] || !(t[1]-t[0]).is_finite())
        { return Err(WindowError::Invalid("invalid or oversized conduction mesh/time grid")); }
        let linear = config.step.linear;
        if !linear.tolerance.is_finite() || linear.tolerance <= 0.0 || linear.tolerance >= 1.0
            || linear.max_iterations == 0 || (config.nonlinear.is_some() && linear.restart == 0)
            || !config.step.energy_tolerance_j.is_finite() || config.step.energy_tolerance_j <= 0.0
        { return Err(WindowError::Invalid("invalid conduction residual/energy policy")); }
        if let Some(policy) = config.nonlinear {
            policy.validate().map_err(|e| refusal(0, "nonlinear policy", e))?;
        }
        // DofMap is the production owner; check foreign boundary indices before
        // calling its indexing constructor. Never silently discard fixed nodes.
        if boundary.dirichlet().iter().any(|&(v,t)| v >= n || !t.is_finite())
            || boundary.dirichlet().windows(2).any(|p| p[0].0 >= p[1].0)
        { return Err(WindowError::Invalid("invalid conduction prescribed-node map")); }
        let dofs = DofMap::new(boundary, n).map_err(|e| refusal(0, "boundary map", e))?;
        let mut owned = zeros(times.len())?; owned.copy_from_slice(times);
        poll(cancelled)?;
        Ok(Self { cx, mesh, dofs, times: owned, config })
    }
    pub fn times(&self) -> &[f64] { &self.times }
    pub fn free_vertices(&self) -> &[usize] { self.dofs.free() }
    pub fn dimension(&self) -> usize { self.dofs.n() }
    pub fn slot_of(&self, vertex: usize) -> Option<usize> {
        if vertex < self.mesh.vertex_count() { self.dofs.slot_of(vertex) } else { None }
    }
    /// Gather a PHYSICAL field, checking its fixed values rather than throwing
    /// them away. Covectors use the private gather path without an affine lift.
    pub fn gather_field(&self, field: &[f64]) -> Result<Vec<f64>, WindowError> {
        if field.len() != self.mesh.vertex_count() || !finite(field)
            || self.dofs.fixed().iter().any(|&v| field[v] != self.dofs.prescribed()[v])
        { return Err(WindowError::Invalid("physical field violates the conduction boundary lift")); }
        Ok(self.dofs.gather(field))
    }
    /// Expand free temperatures with the exact prescribed boundary values.
    pub fn expand_field(&self, free: &[f64]) -> Result<Vec<f64>, WindowError> {
        if free.len() != self.dimension() || !finite(free) {
            return Err(WindowError::Invalid("invalid free temperature field"));
        }
        Ok(self.dofs.scatter(free))
    }
    fn admit_problem(&self, problem: ConductionProblem<'_>) -> Result<(), WindowError> {
        if !std::ptr::eq(problem.mesh, self.mesh)
            || problem.boundary.dirichlet().len() != self.dofs.fixed().len()
            || problem.boundary.dirichlet().iter().zip(self.dofs.fixed())
                .any(|(&(v,t), &expected)| v != expected || t != self.dofs.prescribed()[expected])
        { return Err(WindowError::Invalid("conduction mesh or prescribed boundary changed")); }
        Ok(())
    }
}

/// Complete energy/residual-checked endpoint with its SAME discrete Jacobian.
/// A single interval tape is live during the variational sweep; no PDE trajectory
/// or state covariance is stored here. Endpoint derivative solves use C/dt+J,
/// not a steady operator or differentiation of Newton/Krylov iterations.
pub struct ConductionInterval<'a, 'cx, M> {
    step: StepLinearization<'a>,
    model: &'a M,
    cx: &'a Cx<'cx>,
    dofs: &'a DofMap,
    endpoint: Vec<f64>,
    end: f64,
    interval: usize,
    parameters: usize,
    substep: Option<ConductionSubstep>,
}
impl<M> ConductionInterval<'_, '_, M> {
    /// The actual forward producer report, including discrete energy closure.
    pub fn primal(&self) -> &fs_conduction::transient::backward_euler::StepSolution {
        self.step.primal()
    }
}
impl<'cx, M: ConductionWindowModel> IntervalScheme<M> for ConductionWindowPolicy<'_, 'cx> {
    type Tape<'a> = ConductionInterval<'a, 'cx, M> where Self: 'a, M: 'a;
    fn dimension(&self, _model: &M) -> usize { self.dimension() }
    fn parameter_count(&self, model: &M) -> usize { model.parameter_count() }
    fn validate(&self, times: &[f64]) -> Result<(), WindowError> {
        if times != self.times { return Err(WindowError::Invalid("window must use the declared conduction clock")); }
        Ok(())
    }
    fn record<'a>(&'a self, model: &'a M, interval: usize, start: f64, end: f64,
        initial: &[f64], cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Self::Tape<'a>, WindowError> {
        poll(cancelled)?;
        if interval >= self.times.len()-1 || start != self.times[interval] || end != self.times[interval+1]
            || model.parameter_count() > self.config.max_parameters
        { return Err(WindowError::Invalid("conduction interval or parameter cap mismatch")); }
        self.record_step(model, interval, start, end, initial, None, cancelled)
    }
}
impl<'cx> ConductionWindowPolicy<'_, 'cx> {
    // Shared one-step producer. The public base policy keeps its old data and
    // arithmetic path; refinement supplies only a different fixed time context.
    #[allow(clippy::too_many_arguments)]
    fn record_step<'a, M: ConductionWindowModel>(&'a self, model: &'a M,
        interval: usize, start: f64, end: f64, initial: &[f64],
        substep: Option<ConductionSubstep>, cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<ConductionInterval<'a, 'cx, M>, WindowError> {
        poll(cancelled)?;
        if !start.is_finite() || !end.is_finite() || end <= start || !(end-start).is_finite()
            || model.parameter_count() > self.config.max_parameters
        { return Err(WindowError::Invalid("invalid conduction numerical step")); }
        let old = self.expand_field(initial)?;
        let problem = match substep {
            Some(time) => model.problem_at(time),
            None => model.problem(interval),
        }.map_err(|e| refusal(interval, "endpoint problem", e))?;
        self.admit_problem(problem)?;
        poll(cancelled)?;
        let step = model.engine().linearize_step(self.cx, problem, model.interfaces(), &old,
            end-start, self.config.step, self.config.nonlinear, &[])
            .map_err(|e| refusal(interval, "conduction forward", e))?;
        let endpoint = self.gather_field(&step.primal().temperature)?;
        poll(cancelled)?;
        Ok(ConductionInterval { step, model, cx: self.cx, dofs: &self.dofs,
            endpoint, end, interval, parameters: model.parameter_count(), substep })
    }
}

impl<M: ConductionWindowModel> IntervalTape for ConductionInterval<'_, '_, M> {
    fn endpoint(&self) -> &[f64] { &self.endpoint }
    fn end_time(&self) -> f64 { self.end }
    fn accepted_steps(&self) -> usize { 1 }
    fn pullback(&self, seed: &[f64], direct: &[f64], cancelled: &mut dyn FnMut() -> bool)
        -> Result<TrajectoryGradient, WindowError>
    {
        poll(cancelled)?;
        if seed.len() != self.dofs.n() || direct.len() != self.parameters || !finite(seed) || !finite(direct)
            || self.model.parameter_count() != self.parameters
        { return Err(WindowError::Invalid("invalid conduction endpoint/parameter cotangent")); }
        // A cotangent lift is LINEAR: using scatter would inject the prescribed
        // temperatures as objective weights and corrupt every gradient.
        let mut weights = zeros(self.step.primal().temperature.len())?;
        for (&vertex, &value) in self.dofs.free().iter().zip(seed) { weights[vertex] = value; }
        let gradient = self.step.pullback(self.cx, &weights, &[], &[])
            .map_err(|e| refusal(self.interval, "conduction adjoint", e))?;
        poll(cancelled)?;
        let history = self.step.previous_temperature_pullback(self.cx, &gradient.nodal_load)
            .map_err(|e| refusal(self.interval, "conduction history", e))?;
        let mut parameters = zeros(self.parameters)?; parameters.fill(f64::NAN);
        match self.substep {
            Some(time) => self.model.parameter_pullback_at(time, self.cx, &self.step,
                &gradient.nodal_load, &mut parameters),
            None => self.model.parameter_pullback(self.interval, self.cx, &self.step,
                &gradient.nodal_load, &mut parameters),
        }.map_err(|e| refusal(self.interval, "conduction parameter derivative", e))?;
        poll(cancelled)?;
        if !finite(&parameters) { return Err(WindowError::NonFinite("conduction parameter derivatives")); }
        for (value, partial) in parameters.iter_mut().zip(direct) { *value += partial; }
        if !finite(&parameters) { return Err(WindowError::NonFinite("conduction parameter accumulation")); }
        let initial = self.dofs.gather(&history);
        poll(cancelled)?;
        Ok(TrajectoryGradient { initial, parameters, replayed_steps: 0, peak_checkpoints: 0 })
    }
}

#[cfg(test)]
mod tests;
