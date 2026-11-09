//! Coupled stiff viscosity, Fourier exchange and caloric self-heating.
//!
//! First-order IMEX Euler advances conservative transport/capillarity explicitly,
//! then jointly solves backward Euler for face momentum and extensive internal
//! energy at the transported areas. Every residual uses the actual original
//! caloric temperature, Carreau rate/temperature and conductivity. An original
//! held-temperature momentum solve supplies only a Newton initial guess.
//!
//! The analytic Jacobian is nonsymmetric: thermal softening and viscous heating
//! are coupled, not a positive-definite joint diffusion operator. Native cyclic
//! positive spring solves precondition its frozen momentum and thermal blocks;
//! the existing bounded Newton/FGMRES owner supplies general globalization.
//!
//! Backward Euler's kinetic damping is reported separately and never heated.
//! Stability beyond the explicit diffusion limits does not resolve fast physical
//! heating: maximum_step_s and temporal refinement remain necessary. Advection,
//! capillarity, source support and native geometry limits remain in force.

use std::cell::RefCell;

use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, NonlinearProblem,
};

use super::super::implicit as momentum;
use super::*;

/// Numerical method actually used by this separate coupled advance entry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThermalJetImplicitMethod {
    ImexEuler,
}

/// Bounded nonlinear work and independent physical component tolerances.
#[derive(Clone, Copy, Debug)]
pub struct ThermalJetImplicitControls {
    pub velocity_scale_m_s: f64,
    pub specific_energy_scale_j_kg: f64,
    pub maximum_newton_iterations: usize,
    /// Whole-call joint Newton source-state preparations, including rejected
    /// trials. Transport and the native seed have their own bounded sweeps.
    pub maximum_residual_evaluations: usize,
    pub velocity_tolerance_m_s: f64,
    pub cell_energy_tolerance_j: f64,
    /// Original fixed-temperature momentum solve used only for the initial guess.
    pub seed_momentum: SlenderJetImplicitControls,
    /// Restart and cycle counts bound O(N * restart) native Krylov workspace.
    pub solver: NewtonKrylovConfig,
}
impl Default for ThermalJetImplicitControls {
    fn default() -> Self {
        Self {
            velocity_scale_m_s: 1.0,
            specific_energy_scale_j_kg: 1.0,
            maximum_newton_iterations: 32,
            maximum_residual_evaluations: 100_000,
            velocity_tolerance_m_s: 1e-10,
            cell_energy_tolerance_j: 1e-10,
            seed_momentum: SlenderJetImplicitControls::default(),
            solver: NewtonKrylovConfig::default(),
        }
    }
}
impl ThermalJetImplicitControls {
    fn validate(self) -> Result<(), ThermalJetImplicitError> {
        self.solver
            .validate()
            .map_err(ThermalJetImplicitError::Solver)?;
        self.seed_momentum.validate()?;
        if !self.velocity_scale_m_s.is_finite()
            || self.velocity_scale_m_s <= 0.0
            || !self.specific_energy_scale_j_kg.is_finite()
            || self.specific_energy_scale_j_kg <= 0.0
            || !(1..=1024).contains(&self.maximum_newton_iterations)
            || self.maximum_residual_evaluations == 0
            || !(1..=128).contains(&self.solver.linear_restart)
            || !(1..=1024).contains(&self.solver.max_linear_cycles)
            || !self.velocity_tolerance_m_s.is_finite()
            || self.velocity_tolerance_m_s < 0.0
            || !self.cell_energy_tolerance_j.is_finite()
            || self.cell_energy_tolerance_j < 0.0
        {
            return Err(ThermalJetImplicitError::Invalid(
                "bounded positive scales/work and finite physical tolerances required",
            ));
        }
        Ok(())
    }
}

#[derive(Debug)]
pub enum ThermalJetImplicitError {
    Invalid(&'static str),
    Thermal(ThermalJetError),
    Solver(NewtonError),
    NotConverged(NewtonReport),
    EvaluationBudget,
    ComponentResidual { velocity_m_s: f64, energy_j: f64 },
}
impl From<ThermalJetError> for ThermalJetImplicitError {
    fn from(value: ThermalJetError) -> Self {
        Self::Thermal(value)
    }
}
impl From<NonlinearJetError> for ThermalJetImplicitError {
    fn from(value: NonlinearJetError) -> Self {
        Self::Thermal(value.into())
    }
}
impl core::fmt::Display for ThermalJetImplicitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit thermoviscous jet: {self:?}")
    }
}
impl std::error::Error for ThermalJetImplicitError {}

/// Actual first-order trajectory and distinct physical/numerical accounts.
#[derive(Clone, Copy, Debug)]
pub struct ThermalJetImplicitAdvance {
    pub method: ThermalJetImplicitMethod,
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub requested_duration_s: f64,
    pub accepted_substeps: usize,
    pub attempted_substeps: usize,
    pub outcome: SlenderJetOutcome,
    pub before: SlenderJetDiagnostics,
    pub after: SlenderJetDiagnostics,
    pub viscous_heat_j: f64,
    pub upwind_dissipation_j: f64,
    /// Sum of .5*rho*pi*dx*M_face*(u_end-u_transport)^2 [J].
    pub implicit_kinetic_damping_j: f64,
    pub internal_energy_change_j: f64,
    pub energy_balance_residual_j: f64,
    pub mass_balance_residual_kg: f64,
    /// Delta(kinetic+surface)+physical heat+upwind loss+BE damping.
    /// Explicit transport/capillary truncation and solve defects remain signed.
    pub mechanical_energy_residual_j: f64,
    pub combined_energy_residual_j: f64,
    pub nonlinear_iterations: usize,
    pub linear_iterations: usize,
    pub residual_evaluations: usize,
    pub maximum_velocity_residual_m_s: f64,
    pub maximum_cell_energy_residual_j: f64,
}

#[derive(Clone, Copy, Default)]
struct CellCoefficients {
    mass: f64,
    face_mass: f64,
    cv: f64,
    temperature: f64,
    stress_rate: f64,
    stress_temperature: f64,
    power_rate: f64,
    power_temperature: f64,
}
#[derive(Clone, Copy, Default)]
struct FaceCoefficients {
    conductance: f64,
    left: f64,
    right: f64,
}

pub(super) struct Workspace {
    seed: momentum::ConstitutiveSolve,
    linear: momentum::CyclicSolve,
    cells: Vec<CellCoefficients>,
    faces: Vec<FaceCoefficients>,
    cache: Vec<f64>,
    cache_valid: bool,
    evaluation: Option<Evaluation>,
    input: Vec<f64>,
    output: Vec<f64>,
}
fn filled<T: Clone>(n: usize, value: T, gate: &CancelGate) -> Result<Vec<T>, NonlinearJetError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(n)
        .map_err(|_| NonlinearJetError::Allocation)?;
    while result.len() < n {
        poll(gate)?;
        result.resize(n.min(result.len() + 256), value.clone());
    }
    Ok(result)
}
impl Workspace {
    fn new(n: usize, gate: &CancelGate) -> Result<Self, NonlinearJetError> {
        let count = n.checked_mul(2).ok_or(NonlinearJetError::Allocation)?;
        Ok(Self {
            seed: momentum::ConstitutiveSolve::new(n, gate)?,
            linear: momentum::CyclicSolve::new(n, gate)?,
            cells: filled(n, CellCoefficients::default(), gate)?,
            faces: filled(n, FaceCoefficients::default(), gate)?,
            cache: zeros(count, gate)?,
            cache_valid: false,
            evaluation: None,
            input: zeros(n, gate)?,
            output: zeros(n, gate)?,
        })
    }
}

struct StageContext<'a, 'm> {
    jet: &'a mut ThermalCarreauJet<'m>,
    work: &'a mut Workspace,
    dt: f64,
    controls: ThermalJetImplicitControls,
    evaluations: &'a mut usize,
    gate: &'a CancelGate,
    failure: Option<ThermalJetImplicitError>,
    trial_refusal: Option<ThermalJetImplicitError>,
}
impl StageContext<'_, '_> {
    /// Source values and derivatives at this exact algebraic trial. A cache hit
    /// avoids repeating original material work for Jv and block preconditioning.
    fn prepare(&mut self, x: &[f64]) -> Result<(), ThermalJetImplicitError> {
        poll(self.gate)?;
        let n = self.jet.candidate.len();
        if x.len() != 2 * n || x.iter().any(|v| !v.is_finite()) {
            return Err(ThermalJetImplicitError::Invalid(
                "finite complete stage increments required",
            ));
        }
        if self.work.cache_valid && self.work.cache == x {
            return Ok(());
        }
        if *self.evaluations >= self.controls.maximum_residual_evaluations {
            return Err(ThermalJetImplicitError::EvaluationBudget);
        }
        *self.evaluations += 1;
        self.work.cache_valid = false;
        let rho = self.jet.flow.parameters.density_kg_m3;
        let dx = self.jet.flow.spacing_m;
        let mass_factor = positive(core::f64::consts::PI * rho * dx)?;
        for i in 0..n {
            poll_at(self.gate, i)?;
            let a = self.jet.flow.candidate.a[i];
            let face_mass = positive(average(a, self.jet.flow.candidate.a[next(i, n)]))?;
            let mass = positive(mass_factor * a)?;
            self.jet.flow.stage.a[i] = a;
            self.jet.flow.stage.q[i] = finite(
                self.jet.flow.candidate.q[i] + face_mass * self.controls.velocity_scale_m_s * x[i],
            )?;
            self.jet.stage[i] = finite(
                self.jet.candidate[i] + mass * self.controls.specific_energy_scale_j_kg * x[n + i],
            )?;
        }
        // Sample original physical stresses/heat without constructing an
        // explicit viscous or conductive CFL that this stage removes.
        let evaluation = super::evaluate(
            &self.jet.flow.stage,
            &self.jet.stage,
            &self.jet.source,
            &mut self.jet.workspace,
            self.jet.flow.parameters,
            self.jet.flow.limits,
            self.jet.flow.spacing_m,
            &mut self.jet.flow.spatial,
            &mut self.jet.flow.rhs,
            &mut self.jet.rhs,
            false,
            true,
            self.gate,
        )?;
        momentum::viscous_rhs(
            &self.jet.flow.stage,
            &self.jet.flow.spatial,
            self.jet.flow.parameters,
            Some(&self.jet.workspace.profile),
            dx,
            &mut self.jet.flow.rhs.q,
            self.gate,
        )?;
        for i in 0..n {
            poll_at(self.gate, i)?;
            let temperature = self.jet.workspace.profile.temperature_k(i);
            let cv = self
                .jet
                .source
                .caloric
                .specific_heat_cv_j_kg_k(temperature)
                .map_err(|source| ThermalJetError::Thermal { cell: i, source })?;
            let strain = finite(
                (self.jet.flow.spatial.u[i] - self.jet.flow.spatial.u[previous(i, n)]) / dx,
            )?;
            let sample = self.jet.workspace.profile.sample(strain, i)?;
            let a = self.jet.flow.stage.a[i];
            let mu = sample.dynamic_viscosity_pa_s();
            let tangent = sample.stress_tangent_pa_s();
            let mu_t = sample.temperature_derivative_pa_s_per_k();
            let mechanical = finite(3.0 * a / rho)?;
            let thermal = finite(3.0 * core::f64::consts::PI * dx * a)?;
            self.work.cells[i] = CellCoefficients {
                mass: positive(mass_factor * a)?,
                face_mass: positive(average(a, self.jet.flow.stage.a[next(i, n)]))?,
                cv: positive(cv)?,
                temperature,
                stress_rate: positive(mechanical * tangent)?,
                stress_temperature: finite(mechanical * mu_t * strain)?,
                power_rate: finite(thermal * strain * (mu + tangent))?,
                power_temperature: finite(thermal * mu_t * strain * strain)?,
            };
            if let Some(law) = self.jet.source.conduction {
                let state = law
                    .state_at_temperature(temperature)
                    .map_err(|source| ThermalJetError::Thermal { cell: i, source })?;
                // Reuse scalar scratch for k_T/k, after native flux evaluation.
                self.work.output[i] =
                    finite(state.conductivity_slope_kg_w_m_k_j * cv / state.conductivity_w_m_k)?;
            }
        }
        for i in 0..n {
            poll_at(self.gate, i)?;
            let r = next(i, n);
            self.work.faces[i] = if self.jet.source.conduction.is_some() {
                let left = positive(self.jet.workspace.conductivity[i] * self.jet.flow.stage.a[i])?;
                let right =
                    positive(self.jet.workspace.conductivity[r] * self.jet.flow.stage.a[r])?;
                let scale = left.max(right);
                let weight = (right / scale) / (left / scale + right / scale);
                let g = self.jet.workspace.conductance[i];
                FaceCoefficients {
                    conductance: g,
                    left: finite(g * weight * self.work.output[i])?,
                    right: finite(g * (1.0 - weight) * self.work.output[r])?,
                }
            } else {
                FaceCoefficients::default()
            };
        }
        let powers = self.jet.flow.spatial.cell_viscous_power_w()?;
        for i in 0..n {
            poll_at(self.gate, i)?;
            let l = previous(i, n);
            let flux = self.work.faces[i].conductance
                * (self.work.cells[i].temperature - self.work.cells[next(i, n)].temperature);
            let previous_flux = self.work.faces[l].conductance
                * (self.work.cells[l].temperature - self.work.cells[i].temperature);
            self.jet.rhs[i] = finite(powers[i] - (flux - previous_flux))?;
        }
        self.work.cache.copy_from_slice(x);
        self.work.cache_valid = true;
        self.work.evaluation = Some(evaluation);
        Ok(())
    }

    fn residual(&mut self, x: &[f64], out: &mut [f64]) -> Result<(), ThermalJetImplicitError> {
        self.prepare(x)?;
        let n = self.work.cells.len();
        for i in 0..n {
            poll_at(self.gate, i)?;
            let c = self.work.cells[i];
            out[i] = finite(
                x[i] - self.dt * self.jet.flow.rhs.q[i]
                    / (c.face_mass * self.controls.velocity_scale_m_s),
            )?;
            out[n + i] = finite(
                x[n + i]
                    - self.dt * self.jet.rhs[i]
                        / (c.mass * self.controls.specific_energy_scale_j_kg),
            )?;
        }
        Ok(())
    }

    fn jacobian(
        &mut self,
        x: &[f64],
        direction: &[f64],
        out: &mut [f64],
    ) -> Result<(), ThermalJetImplicitError> {
        self.prepare(x)?;
        let n = self.work.cells.len();
        let dx = self.jet.flow.spacing_m;
        let dv = |i| self.controls.velocity_scale_m_s * direction[i];
        let dt =
            |i| self.controls.specific_energy_scale_j_kg * direction[n + i] / self.work.cells[i].cv;
        let strain = |i| (dv(i) - dv(previous(i, n))) / dx;
        let stress = |i| {
            self.work.cells[i].stress_rate * strain(i)
                + self.work.cells[i].stress_temperature * dt(i)
        };
        let heat_flux = |i| {
            let r = next(i, n);
            let face = self.work.faces[i];
            face.conductance * (dt(i) - dt(r))
                + (face.left * dt(i) + face.right * dt(r))
                    * (self.work.cells[i].temperature - self.work.cells[r].temperature)
        };
        for i in 0..n {
            poll_at(self.gate, i)?;
            let c = self.work.cells[i];
            let gq = (stress(next(i, n)) - stress(i)) / dx;
            let ge = c.power_rate * strain(i) + c.power_temperature * dt(i)
                - (heat_flux(i) - heat_flux(previous(i, n)));
            out[i] = finite(
                direction[i] - self.dt * gq / (c.face_mass * self.controls.velocity_scale_m_s),
            )?;
            out[n + i] = finite(
                direction[n + i]
                    - self.dt * ge / (c.mass * self.controls.specific_energy_scale_j_kg),
            )?;
        }
        Ok(())
    }

    fn precondition(
        &mut self,
        x: &[f64],
        residual: &[f64],
        out: &mut [f64],
    ) -> Result<(), ThermalJetImplicitError> {
        self.prepare(x)?;
        let n = self.work.cells.len();
        let dx = self.jet.flow.spacing_m;
        for i in 0..n {
            poll_at(self.gate, i)?;
            self.work.input[i] = finite(
                self.work.cells[i].face_mass * self.controls.velocity_scale_m_s * residual[i],
            )?;
        }
        let cells = &self.work.cells;
        self.work.linear.solve_weighted(
            &|i| cells[i].face_mass,
            &|i| self.dt * cells[i].stress_rate / dx / dx,
            &self.work.input,
            &mut self.work.output,
            self.gate,
        )?;
        for i in 0..n {
            poll_at(self.gate, i)?;
            out[i] = finite(self.work.output[i] / self.controls.velocity_scale_m_s)?;
            self.work.input[i] =
                finite(cells[i].mass * self.controls.specific_energy_scale_j_kg * residual[n + i])?;
        }
        let faces = &self.work.faces;
        self.work.linear.solve_weighted(
            &|i| cells[i].mass * cells[i].cv,
            &|i| self.dt * faces[previous(i, n)].conductance,
            &self.work.input,
            &mut self.work.output,
            self.gate,
        )?;
        for i in 0..n {
            poll_at(self.gate, i)?;
            out[n + i] = finite(
                self.work.output[i] * cells[i].cv / self.controls.specific_energy_scale_j_kg,
            )?;
        }
        Ok(())
    }

    fn components(&self) -> Result<[f64; 2], ThermalJetImplicitError> {
        let mut velocity = 0.0_f64;
        let mut energy = 0.0_f64;
        for i in 0..self.work.cells.len() {
            poll_at(self.gate, i)?;
            velocity = velocity.max(
                finite(
                    ((self.jet.flow.stage.q[i] - self.jet.flow.candidate.q[i])
                        - self.dt * self.jet.flow.rhs.q[i])
                        / self.work.cells[i].face_mass,
                )?
                .abs(),
            );
            energy = energy.max(
                finite((self.jet.stage[i] - self.jet.candidate[i]) - self.dt * self.jet.rhs[i])?
                    .abs(),
            );
        }
        if velocity > self.controls.velocity_tolerance_m_s
            || energy > self.controls.cell_energy_tolerance_j
        {
            return Err(ThermalJetImplicitError::ComponentResidual {
                velocity_m_s: velocity,
                energy_j: energy,
            });
        }
        Ok([velocity, energy])
    }
}

struct Problem<'a, 'm> {
    context: RefCell<StageContext<'a, 'm>>,
}
impl NonlinearProblem for Problem<'_, '_> {
    fn dimension(&self) -> usize {
        2 * self.context.borrow().work.cells.len()
    }
    fn residual(&self, x: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        let mut c = self.context.borrow_mut();
        if c.failure.is_none() {
            match c.residual(x, out) {
                Ok(()) => c.trial_refusal = None,
                Err(error) => {
                    // Only algebraic line-search candidates may backtrack
                    // within original source support. Accepted-point Jv,
                    // budgets and cancellation retain terminal typed failures.
                    if matches!(
                        &error,
                        ThermalJetImplicitError::Thermal(
                            ThermalJetError::Thermal { .. }
                                | ThermalJetError::Flow(
                                    NonlinearJetError::Rheology { .. }
                                        | NonlinearJetError::NonRepresentable
                                )
                        )
                    ) {
                        c.trial_refusal = Some(error);
                    } else {
                        c.failure = Some(error);
                    }
                    out.fill(f64::NAN);
                }
            }
        }
    }
    fn jacobian_apply(&self, x: &[f64], direction: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        let mut c = self.context.borrow_mut();
        if c.failure.is_none() {
            if let Err(error) = c.jacobian(x, direction, out) {
                c.failure = Some(error);
                out.fill(f64::NAN);
            }
        }
    }
    fn preconditioner_apply(
        &self,
        x: &[f64],
        _: usize,
        _: usize,
        residual: &[f64],
        out: &mut [f64],
    ) {
        out.fill(f64::NAN);
        let mut c = self.context.borrow_mut();
        if c.failure.is_none() {
            if let Err(error) = c.precondition(x, residual, out) {
                c.failure = Some(error);
                out.fill(f64::NAN);
            }
        }
    }
}

struct StageReport {
    evaluation: Evaluation,
    nonlinear: usize,
    linear: usize,
    velocity_residual: f64,
    energy_residual: f64,
}
fn solve_stage(
    jet: &mut ThermalCarreauJet<'_>,
    work: &mut Workspace,
    dt: f64,
    controls: ThermalJetImplicitControls,
    evaluations: &mut usize,
    gate: &CancelGate,
) -> Result<StageReport, ThermalJetImplicitError> {
    let n = jet.candidate.len();
    let mut initial = zeros(n.checked_mul(2).ok_or(NonlinearJetError::Allocation)?, gate)?;
    // Original fixed-T solve is a supported stiff initial guess, not the final
    // physics. Every subsequent Newton residual re-evaluates coupled T and mu.
    for i in 0..n {
        poll_at(gate, i)?;
        jet.flow.stage.a[i] = jet.flow.candidate.a[i];
        jet.flow.stage.q[i] = jet.flow.candidate.q[i];
    }
    if jet.workspace.profile.is_uniform_newtonian() {
        let nu = positive(
            jet.workspace.profile.coefficients(0.0, 0)?[0] / jet.flow.parameters.density_kg_m3,
        )?;
        work.linear
            .apply(&mut jet.flow.stage, nu, jet.flow.spacing_m, dt, gate)?;
    } else {
        work.seed.apply(
            &mut work.linear,
            &mut jet.flow.stage,
            &jet.workspace.profile,
            jet.flow.spacing_m,
            dt,
            controls.seed_momentum,
            gate,
        )?;
    }
    for i in 0..n {
        poll_at(gate, i)?;
        let mass = average(jet.flow.stage.a[i], jet.flow.stage.a[next(i, n)]);
        initial[i] = finite(
            (jet.flow.stage.q[i] - jet.flow.candidate.q[i]) / mass / controls.velocity_scale_m_s,
        )?;
    }
    work.cache_valid = false;
    let problem = Problem {
        context: RefCell::new(StageContext {
            jet,
            work,
            dt,
            controls,
            evaluations,
            gate,
            failure: None,
            trial_refusal: None,
        }),
    };
    let initialized = NewtonKrylovState::new(&problem, initial, controls.solver);
    if let Some(error) = problem.context.borrow_mut().failure.take() {
        return Err(error);
    }
    let mut solver = match initialized {
        Ok(solver) => solver,
        Err(error) => {
            if let Some(source) = problem.context.borrow_mut().trial_refusal.take() {
                return Err(source);
            }
            return Err(ThermalJetImplicitError::Solver(error));
        }
    };
    let result = solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
        problem.context.borrow().failure.is_some() || gate.is_requested()
    });
    if let Some(error) = problem.context.borrow_mut().failure.take() {
        return Err(error);
    }
    poll(gate)?;
    if !result.report.converged {
        if let Some(source) = problem.context.borrow_mut().trial_refusal.take() {
            return Err(source);
        }
        return Err(ThermalJetImplicitError::NotConverged(result.report));
    }
    let mut context = problem.context.borrow_mut();
    context.prepare(&solver.x)?;
    let [velocity_residual, energy_residual] = context.components()?;
    Ok(StageReport {
        evaluation: context
            .work
            .evaluation
            .ok_or(ThermalJetImplicitError::Invalid(
                "missing admitted endpoint",
            ))?,
        nonlinear: result.report.iterations,
        linear: result
            .report
            .history
            .iter()
            .map(|step| step.linear_iterations)
            .sum(),
        velocity_residual,
        energy_residual,
    })
}

impl ThermalCarreauJet<'_> {
    fn evaluate_transport(
        &mut self,
        stage: bool,
        gate: &CancelGate,
    ) -> Result<Evaluation, ThermalJetImplicitError> {
        let fields = if stage {
            &self.flow.stage
        } else {
            &self.flow.candidate
        };
        let energy = if stage { &self.stage } else { &self.candidate };
        Ok(super::evaluate(
            fields,
            energy,
            &self.source,
            &mut self.workspace,
            self.flow.parameters,
            self.flow.limits,
            self.flow.spacing_m,
            &mut self.flow.spatial,
            &mut self.flow.rhs,
            &mut self.rhs,
            true,
            false,
            gate,
        )?)
    }

    /// Advance with first-order IMEX Euler and fully coupled BE viscosity/heat.
    ///
    /// No explicit viscosity or conduction timestep restriction is imposed.
    /// The original transport/capillary CFL, actual-stage source limits and
    /// caller maximum step remain. Newton convergence does not certify temporal
    /// accuracy. Whole-call refusal/cancellation preserves all accepted state.
    /// Native Krylov workspace is O(N*restart); source/preconditioner arrays are
    /// prepared once and reused. This entry does not select the separate native
    /// held-temperature IMEX-SSP3 policy.
    pub fn advance_implicit(
        &mut self,
        duration_s: f64,
        budget: ThermalJetStepBudget,
        controls: ThermalJetImplicitControls,
        gate: &CancelGate,
    ) -> Result<ThermalJetImplicitAdvance, ThermalJetImplicitError> {
        poll(gate)?;
        controls.validate()?;
        let b = budget.flow;
        if !duration_s.is_finite()
            || duration_s < 0.0
            || !b.maximum_step_s.is_finite()
            || b.maximum_step_s <= 0.0
            || !b.cfl.is_finite()
            || b.cfl <= 0.0
            || b.cfl > 0.25
            || b.maximum_substeps == 0
            || !budget.energy_tolerance_j.is_finite()
            || budget.energy_tolerance_j < 0.0
        {
            return Err(ThermalJetImplicitError::Invalid(
                "duration or coupled step budget",
            ));
        }
        let mut work = match self.implicit_workspace.take() {
            Some(work) => work,
            None => Workspace::new(self.accepted.len(), gate)?,
        };
        let result = self.advance_implicit_inner(duration_s, budget, controls, &mut work, gate);
        self.implicit_workspace = Some(work);
        result
    }

    fn advance_implicit_inner(
        &mut self,
        duration_s: f64,
        budget: ThermalJetStepBudget,
        controls: ThermalJetImplicitControls,
        work: &mut Workspace,
        gate: &CancelGate,
    ) -> Result<ThermalJetImplicitAdvance, ThermalJetImplicitError> {
        let target = self.flow.time_s + duration_s;
        if !target.is_finite() || (duration_s > 0.0 && target <= self.flow.time_s) {
            return Err(NonlinearJetError::NonRepresentable.into());
        }
        let n = self.accepted.len();
        for i in 0..n {
            poll_at(gate, i)?;
            self.flow.candidate.a[i] = self.flow.accepted.a[i];
            self.flow.candidate.q[i] = self.flow.accepted.q[i];
            self.candidate[i] = self.accepted[i];
        }
        self.flow.spatial.reset_viscous_heat(gate)?;
        let before = self.flow.diagnostics;
        let b = budget.flow;
        let mut base = self.evaluate_transport(false, gate)?;
        let (mut elapsed, mut cap) = (0.0, b.maximum_step_s);
        let (mut accepted, mut attempted) = (0, 0);
        let (mut heat, mut upwind, mut damping) = (0.0, 0.0, 0.0);
        let (mut evaluations, mut nonlinear, mut linear) = (0, 0, 0);
        let (mut maximum_velocity, mut maximum_energy) = (0.0_f64, 0.0_f64);
        let mut outcome = SlenderJetOutcome::Completed;
        let mass_factor = positive(
            core::f64::consts::PI * self.flow.parameters.density_kg_m3 * self.flow.spacing_m,
        )?;
        while elapsed < duration_s {
            poll(gate)?;
            if attempted >= b.maximum_substeps {
                return Err(NonlinearJetError::SubstepBudget.into());
            }
            let remaining = duration_s - elapsed;
            let dt = positive(remaining.min(cap).min(stable_step(base.rate_per_s, b.cfl)))?;
            let next_elapsed = if dt == remaining {
                duration_s
            } else {
                elapsed + dt
            };
            if next_elapsed <= elapsed
                || self.flow.time_s + next_elapsed <= self.flow.time_s + elapsed
            {
                return Err(NonlinearJetError::NonRepresentable.into());
            }
            attempted += 1;
            for i in 0..n {
                poll_at(gate, i)?;
                self.flow.stage.a[i] = self.flow.candidate.a[i] + dt * self.flow.rhs.a[i];
                self.flow.stage.q[i] = self.flow.candidate.q[i] + dt * self.flow.rhs.q[i];
                self.stage[i] = finite(self.candidate[i] + dt * self.rhs[i])?;
            }
            let transport = match self.evaluate_transport(true, gate) {
                Err(ThermalJetImplicitError::Thermal(ThermalJetError::Flow(
                    NonlinearJetError::InitialDomain(violation),
                ))) => {
                    outcome = SlenderJetOutcome::DomainLimit(violation);
                    break;
                }
                other => other?,
            };
            if dt > stable_step(transport.rate_per_s, b.cfl) {
                cap = positive(0.5 * dt)?;
                base = self.evaluate_transport(false, gate)?;
                continue;
            }
            let numerical = finite(dt * base.upwind_power_w)?;
            // Transported mass, q and energy are the immutable source of this
            // nonlinear stage. Any subsequent Err rolls back the whole call.
            core::mem::swap(&mut self.flow.candidate, &mut self.flow.stage);
            core::mem::swap(&mut self.candidate, &mut self.stage);
            let solved = solve_stage(self, work, dt, controls, &mut evaluations, gate)?;
            let mut step_damping = Sum::default();
            for i in 0..n {
                poll_at(gate, i)?;
                let mass = average(self.flow.stage.a[i], self.flow.stage.a[next(i, n)]);
                let velocity_change =
                    finite((self.flow.stage.q[i] - self.flow.candidate.q[i]) / mass)?;
                step_damping.add(0.5 * mass_factor * mass * velocity_change * velocity_change);
            }
            self.flow.spatial.start_viscous_heat(dt, gate)?;
            self.flow.spatial.accept_viscous_heat_step(gate)?;
            heat = finite(heat + dt * solved.evaluation.viscous_power_w)?;
            upwind = finite(upwind + numerical)?;
            damping = finite(damping + step_damping.total())?;
            nonlinear += solved.nonlinear;
            linear += solved.linear;
            maximum_velocity = maximum_velocity.max(solved.velocity_residual);
            maximum_energy = maximum_energy.max(solved.energy_residual);
            core::mem::swap(&mut self.flow.candidate, &mut self.flow.stage);
            core::mem::swap(&mut self.candidate, &mut self.stage);
            base = self.evaluate_transport(false, gate)?;
            elapsed = next_elapsed;
            accepted += 1;
            cap = b.maximum_step_s;
        }
        if matches!(outcome, SlenderJetOutcome::DomainLimit(_)) {
            base = self.evaluate_transport(false, gate)?;
        }
        let mut energy_change = Sum::default();
        let mut mass_change = Sum::default();
        for i in 0..n {
            poll_at(gate, i)?;
            energy_change.add(self.candidate[i] - self.accepted[i]);
            mass_change.add(mass_factor * (self.flow.candidate.a[i] - self.flow.accepted.a[i]));
        }
        let internal = finite(energy_change.total())?;
        let residual = finite(internal - heat)?;
        if residual.abs() > budget.energy_tolerance_j {
            return Err(ThermalJetError::EnergyBalance {
                residual_j: residual,
                tolerance_j: budget.energy_tolerance_j,
            }
            .into());
        }
        let after = base.diagnostics;
        let mechanical = finite(
            (after.kinetic_energy_j - before.kinetic_energy_j)
                + (after.surface_energy_j - before.surface_energy_j)
                + heat
                + upwind
                + damping,
        )?;
        let end_time = if elapsed == duration_s {
            target
        } else {
            self.flow.time_s + elapsed
        };
        let report = ThermalJetImplicitAdvance {
            method: ThermalJetImplicitMethod::ImexEuler,
            start_time_s: self.flow.time_s,
            end_time_s: end_time,
            requested_duration_s: duration_s,
            accepted_substeps: accepted,
            attempted_substeps: attempted,
            outcome,
            before,
            after,
            viscous_heat_j: heat,
            upwind_dissipation_j: upwind,
            implicit_kinetic_damping_j: damping,
            internal_energy_change_j: internal,
            energy_balance_residual_j: residual,
            mass_balance_residual_kg: finite(mass_change.total())?,
            mechanical_energy_residual_j: mechanical,
            combined_energy_residual_j: finite(mechanical + residual)?,
            nonlinear_iterations: nonlinear,
            linear_iterations: linear,
            residual_evaluations: evaluations,
            maximum_velocity_residual_m_s: maximum_velocity,
            maximum_cell_energy_residual_j: maximum_energy,
        };
        let partition = self.flow.spatial.viscous_heat_partition(heat, gate)?;
        let reference = positive(
            self.workspace.profile.coefficients(0.0, 0)?[0] / self.flow.parameters.density_kg_m3,
        )?;
        poll(gate)?;
        self.flow.spatial.publish_viscous_heat(heat, partition);
        core::mem::swap(&mut self.flow.accepted, &mut self.flow.candidate);
        core::mem::swap(&mut self.accepted, &mut self.candidate);
        if let Some(profile) = &mut self.flow.carreau {
            core::mem::swap(profile, &mut self.workspace.profile);
        }
        self.flow.parameters.kinematic_viscosity_m2_s = reference;
        self.flow.time_s = end_time;
        self.flow.diagnostics = after;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
