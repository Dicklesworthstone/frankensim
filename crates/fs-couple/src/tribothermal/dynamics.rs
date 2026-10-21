//! Finite-body inertia closes the original dry-friction/enthalpy feedback.
//!
//! Two rigid translating bodies slide along the FIRST axis of the contact's
//! unchanged tangent frame. Their inertial masses are the exact receiving
//! thermal cells' masses; those cells represent each complete uniform body.
//! A held normal force maintains the fixed nominal patch without normal work.
//! External tangent forces, rather than prescribed slip, now drive motion.
//!
//! The original fs-tribo return map owns sticking, stopping and signed sliding
//! through `du = du_free - C*traction`, with C=dt²*(1/mA+1/mB). One signed
//! relative-speed unknown makes its source query agree with that solved slip.
//! Every Newton--Krylov trial restarts the original contact/thermal candidate.
//! Both endpoint temperatures, speed and history therefore supply the same
//! material friction query that moves the bodies. Centre-of-mass momentum
//! determines their common motion. Backward-Euler kinetic loss remains a
//! separate numerical account: only original physical friction heats bodies.
//!
//! This retains the native endpoint static-first convention, not an exact
//! within-step stopping-time calculation. A stopping BE step may have zero
//! endpoint slip and finite numerical kinetic loss; refine the interval to
//! resolve physical stopping heat. Static probes require source support at
//! zero speed. No velocity floor, normal geometry update, rotation or moving
//! thermal mesh is introduced.

use super::*;
use fs_thermal::enthalpy::DepositedCell;
use fs_tribo::FrictionRegime;

/// Constant applied forces in the existing first tangent direction, in A/B order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MotionDrive {
    pub normal_load_n: f64,
    pub body_forces_n: [f64; 2],
}

/// Independent momentum and whole-system gates, outside the original heat solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MotionSolve {
    pub sliding: SlidingSolve,
    /// Each trial invokes one complete original contact/thermal candidate solve.
    pub maximum_motion_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Finite-difference perturbation of relative speed [m/s].
    pub difference_step_mps: f64,
    /// Each body's independently reconstructed momentum residual [N s].
    pub impulse_tolerance_n_s: f64,
    /// Includes kinetic storage and the explicit BE numerical-loss account [J].
    pub total_energy_tolerance_j: f64,
    /// Residual is reduced mass times source-minus-solved speed [N s].
    pub newton: NewtonKrylovConfig,
}

#[derive(Clone)]
pub struct SlidingPair<'a> {
    contact: MaterialSlidingContact<'a>,
    inventory: [DepositedCell; 2],
    velocities_mps: [f64; 2],
}

pub struct MotionReport<S> {
    /// Original interface account. Its `drive_work_j` is now the INTERNAL
    /// transfer out of body motion, not an additional exterior supply.
    pub sliding: SlidingReport<S>,
    pub nonlinear: NewtonReport,
    pub motion_trials: usize,
    /// Total original friction evaluations, including all discarded trials.
    pub friction_trials: usize,
    pub body_displacements_m: [f64; 2],
    pub body_momentum_residual_n_s: [f64; 2],
    /// Source-query motion versus the original active set's solved motion.
    pub source_motion_residual_n_s: f64,
    pub kinetic_energy_change_j: f64,
    /// Nonnegative BE kinetic defect, never added to the thermal owner.
    pub numerical_dissipation_j: f64,
    pub external_work_j: f64,
    /// Kinetic + caloric + interface storage + outside friction share + BE loss
    /// minus applied-force work and the original external thermal supply.
    pub total_energy_residual_j: f64,
}

pub struct MotionCandidate<'a, 't, S> {
    thermal: DepositionSimulation<'t>,
    pair: SlidingPair<'a>,
    report: MotionReport<S>,
}
impl<'a, 't, S> MotionCandidate<'a, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub const fn pair(&self) -> &SlidingPair<'a> { &self.pair }
    pub const fn report(&self) -> &MotionReport<S> { &self.report }
    pub fn into_parts(self) -> (DepositionSimulation<'t>, SlidingPair<'a>, MotionReport<S>) {
        (self.thermal, self.pair, self.report)
    }
}

impl<'a> SlidingPair<'a> {
    /// Bind existing finite caloric masses, retaining original contact history.
    /// Velocities are signed components along the FIRST original tangent axis.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        contact: MaterialSlidingContact<'a>,
        velocities_mps: [f64; 2],
    ) -> Result<Self, SlidingError> {
        number(velocities_mps[0] - velocities_mps[1])?;
        if contact.time_s() != thermal.time_s()
            || velocities_mps.iter().any(|v| !v.is_finite()) {
            return Err(SlidingError::Invalid("synchronized finite bodies and finite tangent velocities required"));
        }
        let mut cells = Vec::with_capacity(2);
        for (id, curve) in contact.data.cells.into_iter().zip(contact.data.caloric) {
            let cell = thermal.cells().get(&id)
                .ok_or(SlidingError::Invalid("inactive moving contact body"))?;
            if !thermal.uses_material_law(cell.material_id, curve) {
                return Err(SlidingError::Invalid("moving bodies must retain their original caloric owners"));
            }
            cells.push(cell.clone());
        }
        let pair = Self { contact, inventory: [cells.remove(0), cells.remove(0)], velocities_mps };
        pair.kinetic_energy_j()?;
        Ok(pair)
    }
    pub const fn contact(&self) -> &MaterialSlidingContact<'a> { &self.contact }
    pub const fn velocities_mps(&self) -> [f64; 2] { self.velocities_mps }
    pub fn masses_kg(&self) -> [f64; 2] { self.inventory.each_ref().map(|c| c.mass_kg) }
    pub fn kinetic_energy_j(&self) -> Result<f64, SlidingError> {
        let mut total = 0.0;
        for (cell, velocity) in self.inventory.iter().zip(self.velocities_mps) {
            total = number(total + number(0.5 * cell.mass_kg * velocity * velocity)?)?;
        }
        Ok(total)
    }

    /// Solve both velocities and the original finite thermal network atomically.
    /// The callback has the same other-source/conductance contract as the
    /// original sliding step. Each nonlinear trial restarts every physical
    /// history; `energy_guess_j` uses the original sorted thermal-cell order.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: MotionDrive,
        controls: MotionSolve,
        relative_speed_guess_mps: f64,
        energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MotionCandidate<'a, 't, S>, SlidingError<E>> {
        let cancel = RefCell::new(cancelled);
        let checkpoint = || {
            if (*cancel.borrow_mut())() { Err(SlidingError::Cancelled) } else { Ok(()) }
        };
        checkpoint()?;
        controls.newton.validate()?;
        let dt = number(target_s - self.contact.time_s())?;
        if self.contact.time_s() != thermal.time_s() || dt <= 0.0
            || controls.maximum_motion_trials == 0 || controls.maximum_newton_iterations == 0
            || !relative_speed_guess_mps.is_finite()
            || drive.body_forces_n.iter().any(|v| !v.is_finite())
            || [controls.impulse_tolerance_n_s, controls.total_energy_tolerance_j]
                .iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err(SlidingError::Invalid("invalid motion clocks, forces, controls or sliding guess"));
        }
        for previous in &self.inventory {
            let cell = thermal.cells().get(&previous.id)
                .ok_or(SlidingError::Invalid("moving thermal body disappeared"))?;
            if cell.body_id != previous.body_id || cell.segment_id != previous.segment_id
                || cell.material_id != previous.material_id || cell.mass_kg != previous.mass_kg
                || cell.specific_h_j_kg != previous.specific_h_j_kg
                || cell.geometry != previous.geometry {
                return Err(SlidingError::Invalid("moving body mass, thermal history, identity and reference geometry must be retained"));
            }
        }
        let mass = self.masses_kg();
        let total_mass = number(mass[0] + mass[1])?;
        let momentum = number(number(mass[0] * self.velocities_mps[0])?
            + number(mass[1] * self.velocities_mps[1])?)?;
        let centre = number(number(momentum + number(dt * number(drive.body_forces_n[0]
            + drive.body_forces_n[1])?)?)? / total_mass)?;
        let reduced_mass = number((mass[0] / total_mass) * mass[1])?;
        let compliance = number(dt * number(dt * number(1.0 / mass[0] + 1.0 / mass[1])?)?)?;
        let free_velocity = [number(self.velocities_mps[0] + number(dt * drive.body_forces_n[0] / mass[0])?)?,
            number(self.velocities_mps[1] + number(dt * drive.body_forces_n[1] / mass[1])?)?];
        let released_increment = number(dt * number(free_velocity[0] - free_velocity[1])?)?;
        if reduced_mass <= 0.0 || compliance <= 0.0 {
            return Err(SlidingError::Invalid("finite positive reduced inertia and native compliance required"));
        }
        let trials = Cell::new(0usize);
        let friction_trials = Cell::new(0usize);
        let evaluate = |speed: f64| {
            checkpoint()?;
            number::<E>(speed)?;
            if trials.get() >= controls.maximum_motion_trials { return Err(SlidingError::MotionBudget); }
            trials.set(trials.get() + 1);
            let increment = number(dt * speed)?;
            if speed != 0.0 && increment == 0.0 {
                return Err(SlidingError::Invalid("source-query tangent motion is unrepresentable"));
            }
            let candidate = self.contact.step_candidate_with_motion(thermal, target_s,
                SlidingDrive { normal_load_n: drive.normal_load_n, relative_increment_m: [released_increment, 0.0] },
                compliance, [increment, 0.0], controls.sliding, energy_guess_j,
                |state, _| resolve(state), || (*cancel.borrow_mut())())?;
            friction_trials.set(friction_trials.get().checked_add(candidate.report().friction_trials)
                .ok_or(SlidingError::MotionBudget)?);
            let traction = candidate.report().friction.traction_n;
            let slip = candidate.report().friction.slip_increment_m;
            if traction[1] != 0.0 || slip[1] != 0.0 || number(traction[0] * slip[0])? < 0.0 {
                return Err(SlidingError::Invalid("native contact must retain the passive first-axis response"));
            }
            let relative = number(slip[0] / dt)?;
            let velocity = [number(centre + (mass[1] / total_mass) * relative)?,
                number(centre - (mass[0] / total_mass) * relative)?];
            if slip[0] != 0.0 && (relative == 0.0 || velocity[0] == velocity[1]) {
                return Err(SlidingError::Invalid("native solved body motion is unrepresentable"));
            }
            let residual = [number(number(mass[0] * number(velocity[0] - self.velocities_mps[0])?)?
                    - number(dt * number(drive.body_forces_n[0] - traction[0])?)?)?,
                number(number(mass[1] * number(velocity[1] - self.velocities_mps[1])?)?
                    - number(dt * number(drive.body_forces_n[1] + traction[0])?)?)?];
            let source_residual = number(reduced_mass * number(speed - relative)?)?;
            Ok((candidate, velocity, residual, source_residual))
        };
        // A native zero-slip result is an exact active set, never a speed floor.
        // Use its exact zero source point rather than converging tiny fictitious
        // slip or losing a tiny increment against an existing slip history.
        let zero = evaluate(0.0)?;
        let initial_speed = if zero.0.report().friction.regime == FrictionRegime::Sticking {
            0.0
        } else { relative_speed_guess_mps };
        drop(zero);
        let problem = FallibleFiniteDifference::new(1, controls.difference_step_mps,
            |x: &[f64]| evaluate(x[0]).map(|(_, _, _, residual)| vec![residual]))
            .map_err(numerical)?;
        let construction = NewtonKrylovState::new(&problem, vec![initial_speed], controls.newton);
        if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
        let mut solver = construction?;
        let result = solver.run_controlled(&problem, controls.maximum_newton_iterations,
            &mut || (*cancel.borrow_mut())() || problem.has_failed());
        if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete { return Err(SlidingError::Interrupted(result.progress)); }
        if !result.report.converged { return Err(SlidingError::NotConverged(result.report)); }
        let mut accepted = evaluate(solver.x[0])?;
        if accepted.0.report().friction.regime == FrictionRegime::Sticking && solver.x[0] != 0.0 {
            accepted = evaluate(0.0)?;
            if accepted.0.report().friction.regime != FrictionRegime::Sticking {
                return Err(SlidingError::Invalid("zero-speed source does not admit the proposed native sticking state"));
            }
        }
        let (candidate, velocities, residual, source_residual) = accepted;
        let material = &candidate.report().material;
        let actual_speed = number(candidate.report().friction.slip_increment_m[0] / dt)?.abs();
        let source_speed_residual = number(reduced_mass * number(material.speed_mps - actual_speed)?)?;
        let source_distance_residual = number(number(reduced_mass / dt)?
            * number(material.accumulated_distance_m - candidate.contact().accumulated_distance_m())?)?;
        if residual.iter().any(|r| r.abs() > controls.impulse_tolerance_n_s)
            || [source_residual, source_speed_residual, source_distance_residual]
                .iter().any(|r| r.abs() > controls.impulse_tolerance_n_s) {
            return Err(SlidingError::Invalid("finite-body impulse residual exceeds its independent gate"));
        }
        let mut displacement = [0.0; 2];
        let mut kinetic = 0.0;
        let mut numerical_loss = 0.0;
        let mut work = 0.0;
        for i in 0..2 {
            displacement[i] = number(dt * velocities[i])?;
            let dv = number(velocities[i] - self.velocities_mps[i])?;
            kinetic = number(kinetic + number(0.5 * mass[i] * dv
                * number(velocities[i] + self.velocities_mps[i])?)?)?;
            numerical_loss = number(numerical_loss + number(0.5 * mass[i] * dv * dv)?)?;
            work = number(work + number(drive.body_forces_n[i] * displacement[i])?)?;
        }
        let (thermal, contact, sliding) = candidate.into_parts();
        let total = number(kinetic + sliding.thermal.stored_energy_change_j + sliding.friction.stored_change_j
            + sliding.partition.other_work_j() + numerical_loss - work - sliding.external_heat_j)?;
        if total.abs() > controls.total_energy_tolerance_j {
            return Err(SlidingError::Invalid("moving-body kinetic, thermal, work and numerical-loss accounts do not close"));
        }
        let inventory = self.inventory.each_ref().map(|old| thermal.cells()[&old.id].clone());
        checkpoint()?;
        Ok(MotionCandidate { thermal,
            pair: Self { contact, inventory, velocities_mps: velocities },
            report: MotionReport { sliding, nonlinear: result.report, motion_trials: trials.get(),
                friction_trials: friction_trials.get(), body_displacements_m: displacement,
                body_momentum_residual_n_s: residual, source_motion_residual_n_s: source_residual,
                kinetic_energy_change_j: kinetic,
                numerical_dissipation_j: numerical_loss, external_work_j: work,
                total_energy_residual_j: total } })
    }
}
