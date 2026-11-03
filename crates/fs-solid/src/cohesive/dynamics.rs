//! Small-strain deformable dynamics with retained cohesive interfaces.
//!
//! The implicit midpoint momentum equation uses the material owner's analytic
//! **path-mean** traction, including its surviving-bond viscosity. Its work is
//! therefore conjugate to the actual straight jump path, including prescribed
//! offsets. The endpoint traction tangent is not the derivative of this mean;
//! the existing fallible directional-difference adapter supplies Newton/FGMRES
//! actions at an explicit displacement scale instead.
//!
//! The original reduced mass must be symmetric positive definite and the
//! elastic operator symmetric positive semidefinite. Both must remain pure and
//! unchanged. The bulk potential is `u^T K u/2` in these reduced coordinates;
//! moving essential constraints or additional prestress potentials require
//! their own mechanics/work owner. Reference geometry, trace maps, local frames,
//! areas and source conditions are inherited from `CohesiveState`. Explicit
//! fixed-jump temperature transfer preserves that geometry and all motion while
//! advancing only the original source-defined cohesive coefficients/history.
//! A broken bond keeps its history, mass and velocity, and still resists normal
//! compression through its original law. Explicit damaged-area friction uses
//! that SAME compression/damage path, without a second bond history. No new
//! contact topology, heat routing, finite-strain kinematics or stability/unique-
//! root claim is made.

use std::cell::{Cell, RefCell};
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::cohesive::{CohesiveEnergy, CohesivePoint};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{LinearOp, NewtonKrylovState, NewtonReport, SolverRunProgress};

use super::temperature::{TemperatureConfig, TemperatureEnergy};
use super::{CohesiveState, CohesiveSystem, Config, Error, dot, finite, positive};

/// Source-defined residual friction on the original damaged interface area.
pub mod friction;
use friction::DamagedAreaFriction;

/// Unregularized sourced static/kinetic disks in the original joint solve.
pub mod coulomb;
use coulomb::CoulombBinding;

/// Explicit numerical, trajectory and energy controls. These are not material
/// coefficients or inferred small-strain validity limits.
#[derive(Debug, Clone, Copy)]
pub struct DynamicsConfig {
    /// Reuse the native Newton/Krylov, force, material-energy and displacement
    /// controls. `max_operator_evaluations` counts complete attempted physical
    /// residuals, including both evaluations of each numerical Jacobian action
    /// and the independent final replay. Each residual applies M and K once;
    /// the final energy audit additionally applies each operator twice.
    pub solve: Config,
    /// Maximum-component displacement perturbation of the numerical Jacobian.
    pub difference_step_m: f64,
    pub maximum_step_s: f64,
    /// Bound on each accepted or trial reduced velocity component.
    pub maximum_velocity_m_s: f64,
    pub maximum_steps: usize,
    /// Full kinetic + bulk + interface first-law gate, in joules.
    pub absolute_energy_tolerance_j: f64,
    /// Scales by the magnitudes of endpoint energies, work and physical losses.
    pub relative_energy_tolerance: f64,
}

#[derive(Debug, Clone)]
pub enum DynamicsError {
    Model(Error),
    Energy { residual_j: f64, tolerance_j: f64 },
}
impl From<Error> for DynamicsError {
    fn from(error: Error) -> Self {
        Self::Model(error)
    }
}
impl core::fmt::Display for DynamicsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cohesive dynamics refused: {self:?}")
    }
}
impl std::error::Error for DynamicsError {}

/// Immutable accepted motion. The single original state still owns all
/// displacement and irreversible material history; adopting it never rebonds.
#[derive(Clone)]
pub struct CohesiveMotion<'m, 'a> {
    state: CohesiveState<'m, 'a>,
    mass: &'m dyn LinearOp,
    velocity_m_s: Vec<f64>,
    time_s: f64,
    steps: usize,
    config: DynamicsConfig,
    friction: Option<Arc<[Option<DamagedAreaFriction>]>>,
    coulomb: Option<CoulombBinding>,
}

/// Actual accepted losses at one original bonded quadrature patch, in joules.
/// Fracture loss retains its material meaning; no heat partition is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DynamicsPatchDissipation {
    pub patch_id: u64,
    pub fracture_j: f64,
    pub viscous_j: f64,
    pub friction_j: f64,
}

/// Complete candidate, published only after momentum, material and total-work
/// checks pass. Work/loss fields are incremental, not cumulative heat sources.
pub struct DynamicsStep<'m, 'a> {
    pub motion: CohesiveMotion<'m, 'a>,
    pub nonlinear: NewtonReport,
    pub maximum_force_residual_n: f64,
    pub residual_evaluations: usize,
    /// Original patch order and local axes; includes viscosity and any explicitly
    /// bound damaged-area friction.
    pub mean_traction_pa: Vec<[f64; 3]>,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub bulk_stored_before_j: f64,
    pub bulk_stored_after_j: f64,
    pub cohesive_stored_before_j: f64,
    pub cohesive_stored_after_j: f64,
    /// Rate-independent cohesive work, actual storage change and fracture loss.
    pub cohesive_energy_j: CohesiveEnergy,
    pub viscous_dissipated_j: f64,
    /// Nonnegative work of optional damaged-area friction, separate from bond
    /// viscosity and fracture. No thermal partition is inferred.
    pub friction_dissipated_j: f64,
    /// Same order as the original system patches; computed in the existing
    /// accepted endpoint replay, without another constitutive trial.
    pub patch_dissipation: Vec<DynamicsPatchDissipation>,
    /// Mean reduced external force dotted with the actual displacement change.
    pub external_work_j: f64,
    /// Work entering the interface through changes in prescribed jump offsets.
    /// This includes prescribed viscous and damaged-area friction work once.
    pub prescribed_jump_work_j: f64,
    /// Work input minus storage changes, fracture, viscosity and friction loss.
    pub energy_residual_j: f64,
    pub energy_tolerance_j: f64,
}

/// Fixed-jump temperature candidate with the original inertia and clock.
/// No kinetic/bulk storage changes or mechanical/viscous work occur here.
pub struct DynamicsTemperatureStep<'m, 'a> {
    /// Temperature-transferred interfaces with all original motion retained.
    pub motion: CohesiveMotion<'m, 'a>,
    /// Integrated interface parameter work, storage change and damage release.
    pub energy: TemperatureEnergy,
}

impl<'m, 'a> CohesiveState<'m, 'a> {
    /// Add actual reduced inertia and velocity to this accepted state. The mass
    /// must use exactly the same free DOFs and ordering as the elastic/jump
    /// operators. Symmetry/definiteness are the operator owner's obligations;
    /// checking dimensions and sampled energies cannot establish those globally.
    pub fn into_dynamics(
        self,
        mass: &'m dyn LinearOp,
        velocity_m_s: Vec<f64>,
        time_s: f64,
        config: DynamicsConfig,
    ) -> Result<CohesiveMotion<'m, 'a>, DynamicsError> {
        validate(&self, mass, config)?;
        if velocity_m_s.len() != self.system.n || !time_s.is_finite() || time_s < 0.0 {
            return Err(
                Error::Invalid("complete velocity and finite nonnegative clock required").into(),
            );
        }
        for (&u, &v) in self.displacements.iter().zip(&velocity_m_s) {
            finite(u)?;
            finite(v)?;
            if u.abs() > config.solve.maximum_displacement_m
                || v.abs() > config.maximum_velocity_m_s
            {
                return Err(Error::Invalid("initial motion applicability limit").into());
            }
        }
        Ok(CohesiveMotion {
            state: self,
            mass,
            velocity_m_s,
            time_s,
            steps: 0,
            config,
            friction: None,
            coulomb: None,
        })
    }
}

impl<'m, 'a> CohesiveMotion<'m, 'a> {
    #[must_use]
    pub fn state(&self) -> &CohesiveState<'m, 'a> {
        &self.state
    }
    #[must_use]
    pub fn velocity_m_s(&self) -> &[f64] {
        &self.velocity_m_s
    }
    #[must_use]
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    #[must_use]
    pub fn steps(&self) -> usize {
        self.steps
    }

    /// Change the original sourced bond temperatures without resetting damage,
    /// velocity, time, trajectory count, mass or mechanical step controls.
    /// The next ordinary step evolves under the resulting interface forces.
    pub fn transfer_temperature<'n, 'b>(
        &self,
        target: &'n CohesiveSystem<'b>,
        temperature_axis: &str,
        config: TemperatureConfig,
        cx: &Cx<'_>,
    ) -> Result<DynamicsTemperatureStep<'n, 'b>, DynamicsError>
    where
        'm: 'n,
    {
        self.transfer_temperature_controlled(target, temperature_axis, config, &mut || {
            cx.checkpoint().is_err()
        })
    }

    /// Predicate form; a late refusal leaves the complete original motion intact.
    pub fn transfer_temperature_controlled<'n, 'b, C: FnMut() -> bool>(
        &self,
        target: &'n CohesiveSystem<'b>,
        temperature_axis: &str,
        config: TemperatureConfig,
        cancelled: &mut C,
    ) -> Result<DynamicsTemperatureStep<'n, 'b>, DynamicsError>
    where
        'm: 'n,
    {
        let transferred = self.state.transfer_temperature_controlled(
            target,
            temperature_axis,
            config,
            cancelled,
        )?;
        let friction =
            friction::temperature_successors(self.friction.as_deref(), target, cancelled)?;
        let coulomb = self
            .coulomb
            .as_ref()
            .map(|binding| {
                let next = binding.temperature_successor(target, cancelled)?;
                next.admit_basis(target, self.config)?;
                Ok::<_, Error>(next)
            })
            .transpose()?;
        let mut velocity_m_s = Vec::with_capacity(self.velocity_m_s.len());
        for chunk in self.velocity_m_s.chunks(256) {
            if cancelled() {
                return Err(Error::Cancelled.into());
            }
            velocity_m_s.extend_from_slice(chunk);
        }
        if cancelled() {
            return Err(Error::Cancelled.into());
        }
        Ok(DynamicsTemperatureStep {
            motion: CohesiveMotion {
                state: transferred.state,
                mass: self.mass,
                velocity_m_s,
                time_s: self.time_s,
                steps: self.steps,
                config: self.config,
                friction,
                coulomb,
            },
            energy: transferred.energy,
        })
    }

    /// Advance one straight displacement/jump path. `mean_force_n` is the
    /// applied mean force for this interval, and offsets are endpoint values.
    /// Resolve actual reversals and nonlinear branch selection by time-step
    /// refinement; a converged midpoint step does not certify their resolution.
    pub fn advance(
        &self,
        mean_force_n: &[f64],
        endpoint_prescribed_jumps_m: &[[f64; 3]],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<DynamicsStep<'m, 'a>, DynamicsError> {
        self.advance_controlled(mean_force_n, endpoint_prescribed_jumps_m, dt_s, &mut || {
            cx.checkpoint().is_err()
        })
    }

    /// Predicate form with identical physical evolution. Failed, cancelled or
    /// exhausted trials leave displacement, velocity, clock and damage intact.
    /// Polls every patch and at most 256 vector entries, around opaque operator
    /// calls and inside Newton/FGMRES. Operator apply latency remains its owner's
    /// responsibility. A refused finite-difference stencil refuses the step.
    pub fn advance_controlled<C: FnMut() -> bool>(
        &self,
        mean_force_n: &[f64],
        endpoint_prescribed_jumps_m: &[[f64; 3]],
        dt_s: f64,
        cancelled: &mut C,
    ) -> Result<DynamicsStep<'m, 'a>, DynamicsError> {
        if cancelled() {
            return Err(Error::Cancelled.into());
        }
        validate(&self.state, self.mass, self.config)?;
        positive(dt_s)?;
        if dt_s > self.config.maximum_step_s {
            return Err(Error::Budget("time increment").into());
        }
        let time_s = finite(self.time_s + dt_s)?;
        if time_s <= self.time_s {
            return Err(Error::Invalid("clock increment is unrepresentable").into());
        }
        let steps = self
            .steps
            .checked_add(1)
            .ok_or(Error::Budget("step count overflow"))?;
        if steps > self.config.maximum_steps {
            return Err(Error::Budget("trajectory steps").into());
        }
        if mean_force_n.len() != self.state.system.n
            || endpoint_prescribed_jumps_m.len() != self.state.points.len()
        {
            return Err(Error::Invalid("load or prescribed-jump dimensions differ").into());
        }
        let coulomb_mu = if let Some(binding) = &self.coulomb {
            binding.admit_basis(self.state.system, self.config)?;
            binding.coefficients(self, endpoint_prescribed_jumps_m, dt_s, cancelled)?
        } else {
            Vec::new()
        };
        let problem = StepProblem {
            motion: self,
            force: mean_force_n,
            offsets: endpoint_prescribed_jumps_m,
            dt_s,
            cancelled: RefCell::new(cancelled),
            cancellation_seen: Cell::new(false),
            evaluations: Cell::new(0),
            coulomb_mu: RefCell::new(coulomb_mu),
        };
        for chunk in mean_force_n.chunks(256) {
            problem.poll()?;
            for &value in chunk {
                finite(value)?;
            }
        }
        for offset in endpoint_prescribed_jumps_m {
            problem.poll()?;
            for &value in offset {
                finite(value)?;
            }
        }
        let mut predictor = Vec::with_capacity(self.state.system.n);
        for (i, (&u, &v)) in self
            .state
            .displacements
            .iter()
            .zip(&self.velocity_m_s)
            .enumerate()
        {
            if i % 256 == 0 {
                problem.poll()?;
            }
            predictor.push(finite(u + dt_s * v)?);
        }
        if let Some(binding) = &self.coulomb {
            binding.append_predictor(&mut predictor, &mut || problem.poll().is_err())?;
        }
        let mut attempts = 0;
        let mut history = Vec::new();
        let (solution, nonlinear) = loop {
            let numerical = FallibleFiniteDifference::new(
                predictor.len(),
                self.config.difference_step_m,
                |x: &[f64]| problem.evaluate(x, false).map(|e| e.residual),
            )
            .map_err(difference_error)?;
            let solver = NewtonKrylovState::new(&numerical, predictor, self.config.solve.newton);
            if let Some(error) = numerical.take_failure() {
                return Err(difference_error(error).into());
            }
            problem.poll()?;
            let mut solver = solver.map_err(|e| Error::Solver(e.to_string()))?;
            let run = solver.run_controlled(
                &numerical,
                self.config.solve.max_newton_iterations - attempts,
                &mut || numerical.has_failed() || problem.poll().is_err(),
            );
            if let Some(error) = numerical.take_failure() {
                return Err(difference_error(error).into());
            }
            problem.poll()?;
            if run.progress != SolverRunProgress::Complete {
                return Err(Error::Solver(format!("{:?}", run.progress)).into());
            }
            let mut report = run.report;
            for item in &mut report.history {
                item.iteration += attempts;
            }
            attempts += report.iterations;
            history.append(&mut report.history);
            report.iterations = attempts;
            if !report.converged {
                report.history = history;
                return Err(Error::NotConverged(report).into());
            }
            // Static capture is tested by the same native solve. A patch with
            // resolved slip changes to its original kinetic coefficient, then
            // the complete coupled system is re-solved. Every pass shares the
            // original Newton and physical residual evaluation allowances.
            let changed = if let Some(binding) = &self.coulomb {
                let endpoint_velocity = problem.velocity(&solver.x[..self.state.system.n])?;
                binding.switch_sliding(
                    self,
                    &solver.x[..self.state.system.n],
                    endpoint_prescribed_jumps_m,
                    dt_s,
                    &endpoint_velocity,
                    &mut problem.coulomb_mu.borrow_mut(),
                    &mut || problem.poll().is_err(),
                )?
            } else {
                false
            };
            if !changed {
                report.history = history;
                break (solver.x, report);
            }
            predictor = solver.x;
        };

        // Re-evaluate original physical operators and original accepted histories;
        // the numerical stencil never supplies accepted forces or state.
        let evaluation = problem.evaluate(&solution, true)?;
        let maximum_force_residual_n = evaluation
            .residual
            .iter()
            .fold(0.0_f64, |a, b| a.max(b.abs()));
        if maximum_force_residual_n > self.config.solve.force_tolerance_n {
            return Err(Error::ForceResidual(maximum_force_residual_n).into());
        }
        let endpoint = evaluation
            .endpoint
            .ok_or(Error::Invalid("missing final physical replay"))?;
        let body_solution = &solution[..self.state.system.n];
        let velocity = problem.velocity(body_solution)?;
        let kinetic_before_j = problem.quadratic(self.mass, &self.velocity_m_s)?;
        let kinetic_after_j = problem.quadratic(self.mass, &velocity)?;
        let bulk_stored_before_j =
            problem.quadratic(self.state.system.elastic, &self.state.displacements)?;
        let bulk_stored_after_j = problem.quadratic(self.state.system.elastic, body_solution)?;
        if kinetic_before_j < 0.0 || kinetic_after_j < 0.0 {
            return Err(
                Error::Invalid("mass operator has negative represented kinetic energy").into(),
            );
        }
        let mut external_work_j = 0.0;
        for (i, ((&force, &u1), &u0)) in mean_force_n
            .iter()
            .zip(body_solution)
            .zip(&self.state.displacements)
            .enumerate()
        {
            if i % 256 == 0 {
                problem.poll()?;
            }
            external_work_j = finite(external_work_j + force * (u1 - u0))?;
        }
        let energy_residual_j = finite(
            external_work_j + endpoint.prescribed_work
                - (kinetic_after_j - kinetic_before_j)
                - (bulk_stored_after_j - bulk_stored_before_j)
                - endpoint.energy.stored_change
                - endpoint.energy.dissipated
                - endpoint.viscous
                - endpoint.friction,
        )?;
        let scale = [
            kinetic_before_j,
            kinetic_after_j,
            bulk_stored_before_j,
            bulk_stored_after_j,
            endpoint.stored_before,
            endpoint.stored_after,
            external_work_j,
            endpoint.prescribed_work,
            endpoint.energy.dissipated,
            endpoint.viscous,
            endpoint.friction,
        ]
        .into_iter()
        .try_fold(0.0, |sum, value| finite(sum + value.abs()))?;
        let energy_tolerance_j = finite(
            self.config.absolute_energy_tolerance_j + self.config.relative_energy_tolerance * scale,
        )?;
        if bulk_stored_before_j < -energy_tolerance_j || bulk_stored_after_j < -energy_tolerance_j {
            return Err(
                Error::Invalid("elastic operator has negative represented stored energy").into(),
            );
        }
        if energy_residual_j.abs() > energy_tolerance_j {
            return Err(DynamicsError::Energy {
                residual_j: energy_residual_j,
                tolerance_j: energy_tolerance_j,
            });
        }
        problem.poll()?;
        let mut displacements = solution;
        displacements.truncate(self.state.system.n);
        Ok(DynamicsStep {
            motion: Self {
                state: CohesiveState {
                    system: self.state.system,
                    displacements,
                    points: endpoint.points,
                },
                mass: self.mass,
                velocity_m_s: velocity,
                time_s,
                steps,
                config: self.config,
                friction: self.friction.clone(),
                coulomb: self
                    .coulomb
                    .as_ref()
                    .map(|binding| binding.accepted(endpoint.coulomb_states)),
            },
            nonlinear,
            maximum_force_residual_n,
            residual_evaluations: problem.evaluations.get(),
            mean_traction_pa: endpoint.means,
            kinetic_before_j,
            kinetic_after_j,
            bulk_stored_before_j,
            bulk_stored_after_j,
            cohesive_stored_before_j: endpoint.stored_before,
            cohesive_stored_after_j: endpoint.stored_after,
            cohesive_energy_j: endpoint.energy,
            viscous_dissipated_j: endpoint.viscous,
            friction_dissipated_j: endpoint.friction,
            patch_dissipation: endpoint.patch_dissipation,
            external_work_j,
            prescribed_jump_work_j: endpoint.prescribed_work,
            energy_residual_j,
            energy_tolerance_j,
        })
    }
}

struct StepProblem<'p, 'm, 'a> {
    motion: &'p CohesiveMotion<'m, 'a>,
    force: &'p [f64],
    offsets: &'p [[f64; 3]],
    dt_s: f64,
    cancelled: RefCell<&'p mut dyn FnMut() -> bool>,
    cancellation_seen: Cell<bool>,
    evaluations: Cell<usize>,
    coulomb_mu: RefCell<Vec<f64>>,
}
struct Evaluation<'a> {
    residual: Vec<f64>,
    endpoint: Option<Endpoint<'a>>,
}
struct Endpoint<'a> {
    points: Vec<CohesivePoint<'a>>,
    means: Vec<[f64; 3]>,
    energy: CohesiveEnergy,
    stored_before: f64,
    stored_after: f64,
    viscous: f64,
    friction: f64,
    prescribed_work: f64,
    patch_dissipation: Vec<DynamicsPatchDissipation>,
    coulomb_states: Vec<coulomb::CoulombState>,
}
impl<'p, 'm, 'a> StepProblem<'p, 'm, 'a> {
    fn poll(&self) -> Result<(), Error> {
        if self.cancellation_seen.get() || (self.cancelled.borrow_mut())() {
            self.cancellation_seen.set(true);
            return Err(Error::Cancelled);
        }
        Ok(())
    }
    fn velocity(&self, x: &[f64]) -> Result<Vec<f64>, Error> {
        if x.len() != self.motion.state.system.n {
            return Err(Error::Invalid("trial dimension"));
        }
        let mut velocity = Vec::with_capacity(x.len());
        for (i, ((&u1, &u0), &v0)) in x
            .iter()
            .zip(&self.motion.state.displacements)
            .zip(&self.motion.velocity_m_s)
            .enumerate()
        {
            if i % 256 == 0 {
                self.poll()?;
            }
            finite(u1)?;
            let v1 = finite(2.0 * ((u1 - u0) / self.dt_s) - v0)?;
            if u1.abs() > self.motion.config.solve.maximum_displacement_m
                || v1.abs() > self.motion.config.maximum_velocity_m_s
            {
                return Err(Error::Invalid("trial motion applicability limit"));
            }
            velocity.push(v1);
        }
        Ok(velocity)
    }
    fn apply(&self, operator: &dyn LinearOp, x: &[f64]) -> Result<Vec<f64>, Error> {
        self.poll()?;
        let mut output = vec![f64::NAN; x.len()];
        operator.apply(x, &mut output);
        self.poll()?;
        for chunk in output.chunks(256) {
            self.poll()?;
            for &v in chunk {
                finite(v)?;
            }
        }
        Ok(output)
    }
    fn quadratic(&self, operator: &dyn LinearOp, x: &[f64]) -> Result<f64, Error> {
        let action = self.apply(operator, x)?;
        let mut energy = 0.0;
        for (i, (&a, &b)) in x.iter().zip(&action).enumerate() {
            if i % 256 == 0 {
                self.poll()?;
            }
            energy = finite(energy + 0.5 * a * b)?;
        }
        Ok(energy)
    }
    fn evaluate(&self, x: &[f64], retain: bool) -> Result<Evaluation<'a>, Error> {
        self.poll()?;
        if self.evaluations.get() >= self.motion.config.solve.max_operator_evaluations {
            return Err(Error::Budget("dynamic residual evaluations"));
        }
        self.evaluations.set(self.evaluations.get() + 1);
        let state = &self.motion.state;
        let n = state.system.n;
        let dimension = self
            .motion
            .coulomb
            .as_ref()
            .map_or(Ok(n), |binding| binding.dimension(n))?;
        if x.len() != dimension {
            return Err(Error::Invalid("augmented trial dimension"));
        }
        let coordinates = &x[n..];
        let x = &x[..n];
        let velocity = self.velocity(x)?;
        let mut acceleration = Vec::with_capacity(x.len());
        let mut midpoint = Vec::with_capacity(x.len());
        let mut delta = if retain {
            Vec::with_capacity(x.len())
        } else {
            Vec::new()
        };
        for (i, (((&u1, &u0), &v1), &v0)) in x
            .iter()
            .zip(&state.displacements)
            .zip(&velocity)
            .zip(&self.motion.velocity_m_s)
            .enumerate()
        {
            if i % 256 == 0 {
                self.poll()?;
            }
            acceleration.push(finite((v1 - v0) / self.dt_s)?);
            midpoint.push(finite(0.5 * u0 + 0.5 * u1)?);
            if retain {
                delta.push(finite(u1 - u0)?);
            }
        }
        let mut residual = self.apply(self.motion.mass, &acceleration)?;
        let elastic = self.apply(state.system.elastic, &midpoint)?;
        for (i, ((value, k), force)) in residual.iter_mut().zip(elastic).zip(self.force).enumerate()
        {
            if i % 256 == 0 {
                self.poll()?;
            }
            *value = finite(*value + k - force)?;
        }
        residual.resize(dimension, 0.0);
        let mut endpoint = retain.then(|| Endpoint {
            points: Vec::with_capacity(state.points.len()),
            means: Vec::with_capacity(state.points.len()),
            energy: CohesiveEnergy {
                work: 0.0,
                stored_change: 0.0,
                dissipated: 0.0,
                defect: 0.0,
            },
            stored_before: 0.0,
            stored_after: 0.0,
            viscous: 0.0,
            friction: 0.0,
            prescribed_work: 0.0,
            patch_dissipation: Vec::with_capacity(state.points.len()),
            coulomb_states: Vec::with_capacity(if self.motion.coulomb.is_some() {
                state.points.len()
            } else {
                0
            }),
        });
        for (index, patch) in state.system.patches.iter().enumerate() {
            self.poll()?;
            let old = state.points[index];
            let jump = patch.jump(x, self.offsets[index])?;
            let friction = self
                .motion
                .friction
                .as_ref()
                .and_then(|models| models[index]);
            // Zero mu keeps the exact original material path and arithmetic.
            let (step, mut mean, pressure) = if friction.is_some_and(DamagedAreaFriction::active)
                || self
                    .motion
                    .coulomb
                    .as_ref()
                    .is_some_and(|binding| binding.needs_pressure(index))
            {
                old.trial_over_time_with_damage_weighted_compression(
                    jump,
                    self.dt_s,
                    self.motion.config.solve.local_energy_tolerance_j_m2,
                )
            } else {
                old.trial_over_time_with_mean_traction(
                    jump,
                    self.dt_s,
                    self.motion.config.solve.local_energy_tolerance_j_m2,
                )
                .map(|(step, mean)| (step, mean, 0.0))
            }
            .map_err(|error| Error::Material {
                patch: patch.id,
                error,
            })?;
            let mut friction_work = if let Some(model) = friction.filter(|model| model.active()) {
                let previous = old.jump_m();
                let slip = [
                    finite(jump[1] - previous[1])?,
                    finite(jump[2] - previous[2])?,
                ];
                let (traction, work) = model.response(pressure, slip, self.dt_s)?;
                mean[1] = finite(mean[1] + traction[0])?;
                mean[2] = finite(mean[2] + traction[1])?;
                work
            } else {
                0.0
            };
            if let Some(binding) = &self.motion.coulomb {
                let previous = old.jump_m();
                let slip = [
                    finite(jump[1] - previous[1])?,
                    finite(jump[2] - previous[2])?,
                ];
                let mut traction = [0.0; 2];
                if let Some(c) = binding.coordinate(index) {
                    traction = binding.traction(index, coordinates)?;
                    let cone = binding.residual(
                        patch.id,
                        traction,
                        slip,
                        pressure,
                        self.coulomb_mu.borrow()[index],
                    )?;
                    for d in 0..2 {
                        mean[d + 1] = finite(mean[d + 1] + traction[d])?;
                        residual[n + c + d] = finite(patch.area_m2 * cone[d])?;
                    }
                }
                if let Some(end) = &mut endpoint {
                    let (loss, next) = binding.audit(
                        self.motion,
                        index,
                        self.offsets[index],
                        self.dt_s,
                        &velocity,
                        step.equilibrium.point,
                        traction,
                        slip,
                        pressure,
                        self.coulomb_mu.borrow()[index],
                    )?;
                    friction_work = loss;
                    end.coulomb_states.push(next);
                }
            }
            patch.scatter(mean, &mut residual)?;
            if let Some(end) = &mut endpoint {
                let next = step.equilibrium.point;
                let free_jump = patch.jump(&delta, [0.0; 3])?;
                let prescribed_delta =
                    std::array::from_fn(|c| jump[c] - old.jump_m()[c] - free_jump[c]);
                end.prescribed_work =
                    finite(end.prescribed_work + patch.area_m2 * dot(mean, prescribed_delta))?;
                end.stored_before =
                    finite(end.stored_before + patch.area_m2 * old.stored_energy_j_m2())?;
                end.stored_after =
                    finite(end.stored_after + patch.area_m2 * next.stored_energy_j_m2())?;
                end.energy.work =
                    finite(end.energy.work + patch.area_m2 * step.equilibrium.energy.work)?;
                end.energy.stored_change = finite(
                    end.energy.stored_change
                        + patch.area_m2 * (next.stored_energy_j_m2() - old.stored_energy_j_m2()),
                )?;
                let losses = DynamicsPatchDissipation {
                    patch_id: patch.id,
                    fracture_j: finite(patch.area_m2 * step.equilibrium.energy.dissipated)?,
                    viscous_j: finite(patch.area_m2 * step.viscous_dissipated_j_m2)?,
                    friction_j: finite(patch.area_m2 * friction_work)?,
                };
                if losses.fracture_j < 0.0 || losses.viscous_j < 0.0 || losses.friction_j < 0.0 {
                    return Err(Error::Invalid("negative physical patch dissipation"));
                }
                end.energy.dissipated = finite(end.energy.dissipated + losses.fracture_j)?;
                end.viscous = finite(end.viscous + losses.viscous_j)?;
                end.friction = finite(end.friction + losses.friction_j)?;
                end.patch_dissipation.push(losses);
                end.points.push(next);
                end.means.push(mean);
            }
        }
        if let Some(end) = &mut endpoint {
            end.energy.defect =
                finite(end.energy.work - end.energy.stored_change - end.energy.dissipated)?;
            if end.energy.defect.abs() > self.motion.config.solve.total_cohesive_energy_tolerance_j
            {
                return Err(Error::CohesiveEnergyDefect(end.energy.defect));
            }
        }
        self.poll()?;
        Ok(Evaluation { residual, endpoint })
    }
}

fn difference_error(error: FiniteDifferenceError<Error>) -> Error {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(s) => Error::Invalid(s),
    }
}
fn validate(
    state: &CohesiveState<'_, '_>,
    mass: &dyn LinearOp,
    config: DynamicsConfig,
) -> Result<(), Error> {
    let system = state.system;
    if mass.n() != system.n || system.elastic.n() != system.n {
        return Err(Error::Invalid(
            "unchanged mass and elastic dimensions required",
        ));
    }
    config
        .solve
        .newton
        .validate()
        .map_err(|e| Error::Solver(e.to_string()))?;
    for value in [
        config.difference_step_m,
        config.maximum_step_s,
        config.maximum_velocity_m_s,
        config.absolute_energy_tolerance_j,
        config.solve.force_tolerance_n,
        config.solve.local_energy_tolerance_j_m2,
        config.solve.total_cohesive_energy_tolerance_j,
        config.solve.maximum_displacement_m,
    ] {
        positive(value)?;
    }
    if !config.relative_energy_tolerance.is_finite() || config.relative_energy_tolerance < 0.0 {
        return Err(Error::Invalid(
            "finite nonnegative relative energy tolerance required",
        ));
    }
    if config.maximum_steps == 0
        || config.solve.max_newton_iterations == 0
        || config.solve.max_operator_evaluations == 0
    {
        return Err(Error::Budget("nonzero dynamic work allowances required"));
    }
    let r = config.solve.newton.linear_restart;
    let entries = r
        .checked_add(1)
        .and_then(|rp| {
            system
                .n
                .checked_mul(2)
                .and_then(|n2| n2.checked_add(r))
                .and_then(|width| rp.checked_mul(width))
        })
        .ok_or(Error::Budget("Krylov basis overflow"))?;
    if entries > system.budget.max_krylov_basis_entries {
        return Err(Error::Budget("Krylov basis"));
    }
    Ok(())
}
