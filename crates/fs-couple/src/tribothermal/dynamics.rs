//! Finite-body inertia closes the original dry-friction/enthalpy feedback.
//!
//! Two rigid translating bodies slide along the FIRST axis of the contact's
//! unchanged tangent frame. Their inertial masses are the exact receiving
//! thermal cells' masses; those cells represent each complete uniform body.
//! A held normal force maintains the fixed nominal patch without normal work.
//! External tangent forces, rather than prescribed slip, now drive motion.
//!
//! One relative-speed unknown is solved with the existing Newton--Krylov owner.
//! Every trial calls the original contact/thermal candidate from the same old
//! states. Both endpoint temperatures, speed and slip history therefore supply
//! the same material friction query that slows the bodies. Centre-of-mass
//! momentum determines their common motion. Backward-Euler kinetic loss remains
//! a separate numerical account: only original physical friction heats bodies.
//!
//! This admits strictly positive A-minus-B sliding throughout the interval.
//! A stopping/reversing trial refuses; resolve the stopping event with a contact
//! active-set owner before continuing. There is no velocity floor, sticking
//! reaction, normal/contact geometry update, rotation or moving thermal mesh.

use super::*;
use fs_thermal::enthalpy::DepositedCell;

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
    /// Residual coordinates are impulse, N s; unknown is relative speed, m/s.
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
        let relative = number(velocities_mps[0] - velocities_mps[1])?;
        if contact.time_s() != thermal.time_s()
            || velocities_mps.iter().any(|v| !v.is_finite()) || relative <= 0.0 {
            return Err(SlidingError::Invalid("synchronized finite bodies and positive relative sliding required"));
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
            || !relative_speed_guess_mps.is_finite() || relative_speed_guess_mps <= 0.0
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
        let trials = Cell::new(0usize);
        let friction_trials = Cell::new(0usize);
        let evaluate = |speed: f64| {
            checkpoint()?;
            if !speed.is_finite() || speed <= 0.0 {
                return Err(SlidingError::Invalid("stopping or reversal requires a contact active-set event"));
            }
            if trials.get() >= controls.maximum_motion_trials { return Err(SlidingError::MotionBudget); }
            trials.set(trials.get() + 1);
            let velocity = [number(centre + (mass[1] / total_mass) * speed)?,
                number(centre - (mass[0] / total_mass) * speed)?];
            let relative = number(velocity[0] - velocity[1])?;
            let increment = number(dt * relative)?;
            if relative <= 0.0 || increment <= 0.0 {
                return Err(SlidingError::Invalid("positive relative body motion is unrepresentable"));
            }
            let candidate = self.contact.step_candidate(thermal, target_s,
                SlidingDrive { normal_load_n: drive.normal_load_n, relative_increment_m: [increment, 0.0] },
                controls.sliding, energy_guess_j, &resolve, || (*cancel.borrow_mut())())?;
            friction_trials.set(friction_trials.get().checked_add(candidate.report().friction_trials)
                .ok_or(SlidingError::MotionBudget)?);
            let traction = candidate.report().friction.traction_n;
            if traction[0] < 0.0 || traction[1] != 0.0 {
                return Err(SlidingError::Invalid("original traction must oppose the admitted positive slide"));
            }
            let residual = [number(number(mass[0] * number(velocity[0] - self.velocities_mps[0])?)?
                    - number(dt * number(drive.body_forces_n[0] - traction[0])?)?)?,
                number(number(mass[1] * number(velocity[1] - self.velocities_mps[1])?)?
                    - number(dt * number(drive.body_forces_n[1] + traction[0])?)?)?];
            Ok((candidate, velocity, residual))
        };
        let problem = FallibleFiniteDifference::new(1, controls.difference_step_mps,
            |x: &[f64]| evaluate(x[0]).map(|(_, _, residual)| vec![residual[0]]))
            .map_err(numerical)?;
        let construction = NewtonKrylovState::new(&problem, vec![relative_speed_guess_mps], controls.newton);
        if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
        let mut solver = construction?;
        let result = solver.run_controlled(&problem, controls.maximum_newton_iterations,
            &mut || (*cancel.borrow_mut())() || problem.has_failed());
        if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete { return Err(SlidingError::Interrupted(result.progress)); }
        if !result.report.converged { return Err(SlidingError::NotConverged(result.report)); }
        let (candidate, velocities, residual) = evaluate(solver.x[0])?;
        if residual.iter().any(|r| r.abs() > controls.impulse_tolerance_n_s) {
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
                body_momentum_residual_n_s: residual, kinetic_energy_change_j: kinetic,
                numerical_dissipation_j: numerical_loss, external_work_j: work,
                total_energy_residual_j: total } })
    }
}
