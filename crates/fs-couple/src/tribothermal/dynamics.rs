//! Finite-body inertia closes the original dry-friction/enthalpy feedback.
//!
//! Two rigid translating bodies slide in both axes of the contact's unchanged
//! tangent frame. Material-bound pairs use the original resolved-body
//! inertias; each receiving thermal cell must contain that complete uniform
//! body's material mass and volume.
//! A held normal force maintains the fixed nominal patch without normal work.
//! External tangent forces, rather than prescribed slip, now drive motion.
//!
//! The original fs-tribo return map owns sticking, stopping and signed sliding
//! through `du = du_free - C*traction`, with C=dt²*(1/mA+1/mB) times the tangent
//! identity. Two signed relative-velocity components make its source query
//! agree with that solved slip; the original material owner queries its norm.
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
//! zero speed. No velocity floor, normal geometry update, rotation, wear or
//! moving thermal mesh is introduced.

/// Work-conjugate collinear motion with physical stopping and reversal events.
pub mod collinear;

use super::*;
use fs_thermal::enthalpy::DepositedCell;
use fs_mbd_material::ResolvedBody;
use std::sync::Arc;
use fs_tribo::FrictionRegime;

/// Constant applied tangent forces, in the existing fixed contact frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MotionDrive {
    pub normal_load_n: f64,
    /// `[[A tangent0, A tangent1], [B tangent0, B tangent1]]` [N].
    pub body_forces_n: [[f64; 2]; 2],
}

/// Independent momentum and whole-system gates, outside the original heat solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MotionSolve {
    pub sliding: SlidingSolve,
    /// Each trial invokes one complete original contact/thermal candidate solve.
    pub maximum_motion_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Finite-difference perturbation of each relative-velocity component [m/s].
    pub difference_step_mps: f64,
    /// Euclidean norm of each body's independently reconstructed momentum residual [N s].
    pub impulse_tolerance_n_s: f64,
    /// Includes kinetic storage and the explicit BE numerical-loss account [J].
    pub total_energy_tolerance_j: f64,
    /// Residual is reduced mass times source-minus-solved velocity [N s].
    pub newton: NewtonKrylovConfig,
}

#[derive(Clone)]
pub struct SlidingPair<'a> {
    contact: MaterialSlidingContact<'a>,
    inventory: [DepositedCell; 2],
    material_bodies: Option<Arc<[ResolvedBody; 2]>>,
    velocities_mps: [[f64; 2]; 2],
}

pub struct MotionReport<S> {
    /// Original interface account. Its `drive_work_j` is now the INTERNAL
    /// transfer out of body motion, not an additional exterior supply.
    pub sliding: SlidingReport<S>,
    pub nonlinear: NewtonReport,
    pub motion_trials: usize,
    /// Total original friction evaluations, including all discarded trials.
    pub friction_trials: usize,
    /// A/B body order, then the two fixed tangent-frame components.
    pub body_displacements_m: [[f64; 2]; 2],
    /// A/B body order, then the two fixed tangent-frame components.
    pub body_momentum_residual_n_s: [[f64; 2]; 2],
    /// Source-query motion versus the original active set's solved motion.
    pub source_motion_residual_n_s: [f64; 2],
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
    /// Bind explicitly authored finite caloric masses, retaining contact history.
    /// Use `from_material_bodies` to consume source-resolved density and geometry.
    /// Velocities are `[[A tangent0, A tangent1], [B tangent0, B tangent1]]`.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        contact: MaterialSlidingContact<'a>,
        velocities_mps: [[f64; 2]; 2],
    ) -> Result<Self, SlidingError> {
        for axis in 0..2 { number(velocities_mps[0][axis] - velocities_mps[1][axis])?; }
        if contact.time_s() != thermal.time_s()
            || velocities_mps.iter().flatten().any(|v| !v.is_finite()) {
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
        let pair = Self { contact, inventory: [cells.remove(0), cells.remove(0)], material_bodies: None, velocities_mps };
        pair.kinetic_energy_j()?;
        Ok(pair)
    }
    /// Bind source-resolved rigid bodies to their ordered, complete caloric cells.
    ///
    /// Each body's regions must use the original surface material card. Total
    /// source mass and volume must match the finite inventory to roundoff
    /// (128 machine epsilons relative); this is not a material-uncertainty allowance.
    /// The original immutable bodies and density receipts survive every step.
    /// Heating preserves reference mass; it does not re-evaluate density at
    /// unchanged geometry or recreate a live body's momentum.
    pub fn from_material_bodies(
        thermal: &DepositionSimulation<'_>,
        contact: MaterialSlidingContact<'a>,
        bodies: [ResolvedBody; 2],
        velocities_mps: [[f64; 2]; 2],
    ) -> Result<Self, SlidingError> {
        let mut pair = Self::new(thermal, contact, velocities_mps)?;
        for (index, body) in bodies.iter().enumerate() {
            let identity = pair.contact.binding().point.materials[index].content_hash();
            if body.regions().iter().any(|r| r.material().card_identity() != identity) {
                return Err(SlidingError::Invalid("moving body must retain its original ordered surface material"));
            }
            let volume = body.regions().iter().try_fold(0.0,
                |sum, region| number::<Infallible>(sum + region.geometry().volume_m3()))?;
            let cell = &pair.inventory[index];
            if !same_inventory_quantity(body.inertia().mass(), cell.mass_kg)
                || !same_inventory_quantity(volume, cell.geometry.volume_m3) {
                return Err(SlidingError::Invalid("sourced body mass and volume must match its complete caloric cell"));
            }
        }
        pair.material_bodies = Some(Arc::new(bodies));
        pair.kinetic_energy_j()?;
        Ok(pair)
    }
    /// Original density, geometry and material receipts, when source-bound.
    pub fn material_bodies(&self) -> Option<&[ResolvedBody; 2]> {
        self.material_bodies.as_deref()
    }
    pub const fn contact(&self) -> &MaterialSlidingContact<'a> { &self.contact }
    /// A/B body order, then both components of the unchanged tangent frame.
    pub const fn velocities_mps(&self) -> [[f64; 2]; 2] { self.velocities_mps }
    pub fn masses_kg(&self) -> [f64; 2] {
        match self.material_bodies.as_deref() {
            Some(bodies) => bodies.each_ref().map(|body| body.inertia().mass()),
            None => self.inventory.each_ref().map(|cell| cell.mass_kg),
        }
    }
    pub fn kinetic_energy_j(&self) -> Result<f64, SlidingError> {
        let mut total = 0.0;
        for (mass, velocity) in self.masses_kg().into_iter().zip(self.velocities_mps) {
            for component in velocity {
                total = number(total + number(0.5 * mass * component * component)?)?;
            }
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
        relative_velocity_guess_mps: [f64; 2],
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
            || relative_velocity_guess_mps.iter().any(|v| !v.is_finite())
            || drive.body_forces_n.iter().flatten().any(|v| !v.is_finite())
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
        let reduced_mass = number((mass[0] / total_mass) * mass[1])?;
        let compliance = number(dt * number(dt * number(1.0 / mass[0] + 1.0 / mass[1])?)?)?;
        let mut centre = [0.0; 2];
        let mut released_increment = [0.0; 2];
        for axis in 0..2 {
            let momentum = number(number(mass[0] * self.velocities_mps[0][axis])?
                + number(mass[1] * self.velocities_mps[1][axis])?)?;
            centre[axis] = number(number(momentum + number(dt * number(drive.body_forces_n[0][axis]
                + drive.body_forces_n[1][axis])?)?)? / total_mass)?;
            let free_velocity = [number(self.velocities_mps[0][axis]
                    + number(dt * drive.body_forces_n[0][axis] / mass[0])?)?,
                number(self.velocities_mps[1][axis]
                    + number(dt * drive.body_forces_n[1][axis] / mass[1])?)?];
            released_increment[axis] = number(dt * number(free_velocity[0] - free_velocity[1])?)?;
        }
        if reduced_mass <= 0.0 || compliance <= 0.0 {
            return Err(SlidingError::Invalid("finite positive reduced inertia and native compliance required"));
        }
        let trials = Cell::new(0usize);
        let friction_trials = Cell::new(0usize);
        let evaluate = |source_velocity: [f64; 2]| {
            checkpoint()?;
            for component in source_velocity { number::<E>(component)?; }
            if trials.get() >= controls.maximum_motion_trials { return Err(SlidingError::MotionBudget); }
            trials.set(trials.get() + 1);
            let mut increment = [0.0; 2];
            for axis in 0..2 {
                increment[axis] = number(dt * source_velocity[axis])?;
                if source_velocity[axis] != 0.0 && increment[axis] == 0.0 {
                    return Err(SlidingError::Invalid("source-query tangent motion is unrepresentable"));
                }
            }
            let candidate = self.contact.step_candidate_with_motion(thermal, target_s,
                SlidingDrive { normal_load_n: drive.normal_load_n, relative_increment_m: released_increment },
                compliance, increment, controls.sliding, energy_guess_j,
                |state, _| resolve(state), || (*cancel.borrow_mut())())?;
            friction_trials.set(friction_trials.get().checked_add(candidate.report().friction_trials)
                .ok_or(SlidingError::MotionBudget)?);
            let traction = candidate.report().friction.traction_n;
            let slip = candidate.report().friction.slip_increment_m;
            let friction_work = number(number(traction[0] * slip[0])?
                + number(traction[1] * slip[1])?)?;
            if friction_work < 0.0 {
                return Err(SlidingError::Invalid("native dry contact must retain passive tangent work"));
            }
            let mut velocity = [[0.0; 2]; 2];
            let mut residual = [[0.0; 2]; 2];
            let mut source_residual = [0.0; 2];
            for axis in 0..2 {
                let relative = number(slip[axis] / dt)?;
                velocity[0][axis] = number(centre[axis] + (mass[1] / total_mass) * relative)?;
                velocity[1][axis] = number(centre[axis] - (mass[0] / total_mass) * relative)?;
                if slip[axis] != 0.0 && (relative == 0.0 || velocity[0][axis] == velocity[1][axis]) {
                    return Err(SlidingError::Invalid("native solved body motion is unrepresentable"));
                }
                residual[0][axis] = number(number(mass[0]
                        * number(velocity[0][axis] - self.velocities_mps[0][axis])?)?
                    - number(dt * number(drive.body_forces_n[0][axis] - traction[axis])?)?)?;
                residual[1][axis] = number(number(mass[1]
                        * number(velocity[1][axis] - self.velocities_mps[1][axis])?)?
                    - number(dt * number(drive.body_forces_n[1][axis] + traction[axis])?)?)?;
                source_residual[axis] = number(reduced_mass * number(source_velocity[axis] - relative)?)?;
            }
            Ok((candidate, velocity, residual, source_residual))
        };
        // A native zero-slip result is an exact active set, never a speed floor.
        // Use its exact zero source point rather than converging tiny fictitious
        // slip or losing a tiny increment against an existing slip history.
        let zero = evaluate([0.0; 2])?;
        let initial_velocity = if zero.0.report().friction.regime == FrictionRegime::Sticking {
            [0.0; 2]
        } else { relative_velocity_guess_mps };
        drop(zero);
        let problem = FallibleFiniteDifference::new(2, controls.difference_step_mps,
            |x: &[f64]| evaluate([x[0], x[1]]).map(|(_, _, _, residual)| residual.to_vec()))
            .map_err(numerical)?;
        let construction = NewtonKrylovState::new(&problem, initial_velocity.to_vec(), controls.newton);
        if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
        let mut solver = construction?;
        let result = solver.run_controlled(&problem, controls.maximum_newton_iterations,
            &mut || (*cancel.borrow_mut())() || problem.has_failed());
        if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete { return Err(SlidingError::Interrupted(result.progress)); }
        if !result.report.converged { return Err(SlidingError::NotConverged(result.report)); }
        let solved_velocity = [solver.x[0], solver.x[1]];
        let mut accepted = evaluate(solved_velocity)?;
        if accepted.0.report().friction.regime == FrictionRegime::Sticking && solved_velocity != [0.0; 2] {
            accepted = evaluate([0.0; 2])?;
            if accepted.0.report().friction.regime != FrictionRegime::Sticking {
                return Err(SlidingError::Invalid("zero-speed source does not admit the proposed native sticking state"));
            }
        }
        let (candidate, velocities, residual, source_residual) = accepted;
        let material = &candidate.report().material;
        let actual_slip = candidate.report().friction.slip_increment_m;
        let actual_speed = number(number(det::hypot(actual_slip[0], actual_slip[1]))? / dt)?;
        let source_speed_residual = number(reduced_mass * number(material.speed_mps - actual_speed)?)?;
        let source_distance_residual = number(number(reduced_mass / dt)?
            * number(material.accumulated_distance_m - candidate.contact().accumulated_distance_m())?)?;
        if residual.iter().any(|r| det::hypot(r[0], r[1]) > controls.impulse_tolerance_n_s)
            || [number(det::hypot(source_residual[0], source_residual[1]))?,
                source_speed_residual, source_distance_residual]
                .iter().any(|r| r.abs() > controls.impulse_tolerance_n_s) {
            return Err(SlidingError::Invalid("finite-body impulse residual exceeds its independent gate"));
        }
        let mut displacement = [[0.0; 2]; 2];
        let mut kinetic = 0.0;
        let mut numerical_loss = 0.0;
        let mut work = 0.0;
        for i in 0..2 {
            for axis in 0..2 {
                displacement[i][axis] = number(dt * velocities[i][axis])?;
                let dv = number(velocities[i][axis] - self.velocities_mps[i][axis])?;
                kinetic = number(kinetic + number(0.5 * mass[i] * dv
                    * number(velocities[i][axis] + self.velocities_mps[i][axis])?)?)?;
                numerical_loss = number(numerical_loss + number(0.5 * mass[i] * dv * dv)?)?;
                work = number(work + number(drive.body_forces_n[i][axis] * displacement[i][axis])?)?;
            }
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
            pair: Self { contact, inventory, material_bodies: self.material_bodies.clone(), velocities_mps: velocities },
            report: MotionReport { sliding, nonlinear: result.report, motion_trials: trials.get(),
                friction_trials: friction_trials.get(), body_displacements_m: displacement,
                body_momentum_residual_n_s: residual, source_motion_residual_n_s: source_residual,
                kinetic_energy_change_j: kinetic,
                numerical_dissipation_j: numerical_loss, external_work_j: work,
                total_energy_residual_j: total } })
    }
}

fn same_inventory_quantity(source: f64, inventory: f64) -> bool {
    (source - inventory).abs() <= 128.0 * f64::EPSILON * source.abs().max(inventory.abs())
}
