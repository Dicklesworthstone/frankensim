//! Inertial two-body normal motion with the original finite-caloric contact.
//!
//! The single outer unknown is endpoint indentation. The original normal owner
//! integrates Hertz/Hunt--Crossley force and physical work over its discrete
//! indentation segment. That SAME mean force supplies opposite body impulses;
//! actual midpoint travel closes kinetic work without a numerical heat source.
//! Each trial restarts accepted mechanical and thermal histories. Temperature,
//! normal stiffness, load and patch geometry respond inside the nested solve.
//!
//! Each caloric cell is one complete translating body. The normal direction and
//! reference geometry are fixed. This is a finite-contact branch, not an impact
//! or flight event locator: shorten intervals that reach separation or leave the
//! original no-tension/rate/source domains. Segment-rate damping and thermal
//! composition need time refinement, especially across velocity reversal.

use super::*;
use fs_thermal::enthalpy::DepositedCell;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalMotionSolve {
    pub contact: NormalThermalSolve,
    pub maximum_motion_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Perturbation of the endpoint indentation [m].
    pub difference_step_m: f64,
    pub displacement_tolerance_m: f64,
    pub momentum_tolerance_n_s: f64,
    pub total_energy_tolerance_j: f64,
    /// One unknown and one midpoint-closure residual, both in meters.
    pub newton: NewtonKrylovConfig,
}

/// Velocities share one fixed normal axis, directed from body A toward B.
/// A-minus-B motion increases indentation; compressive force acts negatively
/// on A and positively on B. Thermal cells retain the complete body masses.
#[derive(Clone)]
pub struct InertialNormalPair<'a, M: NormalModulus = AffineNormalModulus> {
    contact: CaloricNormalContact<'a, M>,
    inventory: [DepositedCell; 2],
    velocities_mps: [f64; 2],
    translations_m: [f64; 2],
}

pub struct NormalMotionReport<S> {
    pub contact: NormalThermalReport<S>,
    pub nonlinear: NewtonReport,
    pub motion_trials: usize,
    /// Includes successful thermal solves discarded by the outer iteration.
    pub contact_trials: usize,
    pub body_displacements_m: [f64; 2],
    pub relative_midpoint_velocity_mps: f64,
    /// Positive compressive impulse acting positively on B, negatively on A.
    pub compressive_impulse_n_s: f64,
    pub body_momentum_residual_n_s: [f64; 2],
    pub indentation_residual_m: f64,
    pub kinetic_energy_change_j: f64,
    pub external_work_j: f64,
    /// Delta K + native normal work - external mechanical work.
    pub mechanical_energy_residual_j: f64,
    /// Delta(K + U_normal + H_caloric) - external work - other heat.
    pub total_energy_residual_j: f64,
}

pub struct NormalMotionCandidate<'a, 't, S, M: NormalModulus = AffineNormalModulus> {
    thermal: DepositionSimulation<'t>,
    pair: InertialNormalPair<'a, M>,
    report: NormalMotionReport<S>,
}
impl<'a, 't, S, M: NormalModulus> NormalMotionCandidate<'a, 't, S, M> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn pair(&self) -> &InertialNormalPair<'a, M> {
        &self.pair
    }
    pub const fn report(&self) -> &NormalMotionReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        InertialNormalPair<'a, M>,
        NormalMotionReport<S>,
    ) {
        (self.thermal, self.pair, self.report)
    }
}
impl<'a, M: NormalModulus> InertialNormalPair<'a, M> {
    /// Attach without resetting loaded contact or caloric history. Translations
    /// are measured from attachment; supplied velocities use the common axis.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        contact: CaloricNormalContact<'a, M>,
        velocities_mps: [f64; 2],
    ) -> Result<Self, NormalThermalError> {
        if thermal.time_s() != contact.time_s()
            || temperatures(thermal, contact.curves, contact.cells)? != contact.temperature_k()
            || velocities_mps.iter().any(|v| !v.is_finite())
        {
            return Err(NormalThermalError::Invalid(
                "synchronized normal/caloric state and finite velocities required",
            ));
        }
        let inventory = contact.cells.map(|id| thermal.cells()[&id].clone());
        whole_bodies(thermal, &inventory)?;
        if inventory.each_ref().map(|cell| cell.mass_kg) != contact.mass_kg
            || inventory
                .iter()
                .any(|cell| !cell.mass_kg.is_finite() || cell.mass_kg <= 0.0)
        {
            return Err(NormalThermalError::Invalid(
                "positive original complete body masses required",
            ));
        }
        let pair = Self {
            contact,
            inventory,
            velocities_mps,
            translations_m: [0.0; 2],
        };
        pair.kinetic_energy_j()?;
        Ok(pair)
    }
    pub const fn contact(&self) -> &CaloricNormalContact<'a, M> {
        &self.contact
    }
    pub const fn velocities_mps(&self) -> [f64; 2] {
        self.velocities_mps
    }
    pub const fn translations_m(&self) -> [f64; 2] {
        self.translations_m
    }
    pub fn masses_kg(&self) -> [f64; 2] {
        self.inventory.each_ref().map(|cell| cell.mass_kg)
    }
    pub fn kinetic_energy_j(&self) -> Result<f64, NormalThermalError> {
        let mut energy = 0.0;
        for (mass, v) in self.masses_kg().into_iter().zip(self.velocities_mps) {
            energy = finite(energy + finite(0.5 * mass * v * v)?)?;
        }
        Ok(energy)
    }

    /// Constant exterior normal forces drive both actual inertias. The contact
    /// and thermal solve is repeated at each proposed endpoint indentation.
    /// `energy_guess_j` retains the original thermal owner's increment ordering;
    /// callback supplies only OTHER heat/conductances at its proposed endpoint.
    /// At most maximum_motion_trials native solves run, each with its explicit
    /// inner budget. Source failure is propagated, never replaced by a trial.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        body_forces_n: [f64; 2],
        controls: NormalMotionSolve,
        indentation_guess_m: f64,
        energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<NormalMotionCandidate<'a, 't, S, M>, NormalThermalError<E>> {
        let cancel = RefCell::new(cancelled);
        let poll = || {
            if (*cancel.borrow_mut())() {
                Err(NormalThermalError::Cancelled)
            } else {
                Ok(())
            }
        };
        poll()?;
        controls
            .newton
            .validate()
            .map_err(NormalThermalError::Newton)?;
        let dt = finite(target_s - self.contact.time_s())?;
        if thermal.time_s() != self.contact.time_s()
            || dt <= 0.0
            || controls.maximum_motion_trials == 0
            || controls.maximum_newton_iterations == 0
            || !indentation_guess_m.is_finite()
            || indentation_guess_m < 0.0
            || body_forces_n.iter().any(|v| !v.is_finite())
            || [
                controls.displacement_tolerance_m,
                controls.momentum_tolerance_n_s,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(NormalThermalError::Invalid(
                "invalid normal-motion state, drive or controls",
            ));
        }
        for old in &self.inventory {
            let now = thermal
                .cells()
                .get(&old.id)
                .ok_or(NormalThermalError::Invalid(
                    "complete normal body disappeared",
                ))?;
            if now.body_id != old.body_id
                || now.segment_id != old.segment_id
                || now.material_id != old.material_id
                || now.mass_kg != old.mass_kg
                || now.geometry != old.geometry
                || now.specific_h_j_kg != old.specific_h_j_kg
            {
                return Err(NormalThermalError::Invalid(
                    "original normal body inventory or accepted thermal state changed",
                ));
            }
        }
        whole_bodies(thermal, &self.inventory)?;
        let masses = self.masses_kg();
        let inverse_mass = finite(masses[0].recip() + masses[1].recip())?;
        let free_relative_midpoint = finite(
            self.velocities_mps[0] - self.velocities_mps[1]
                + finite(
                    0.5 * dt * finite(body_forces_n[0] / masses[0] - body_forces_n[1] / masses[1])?,
                )?,
        )?;
        let old_indentation = self.contact.endpoint().approach_m;
        let free_indentation = finite(old_indentation + finite(dt * free_relative_midpoint)?)?;
        if free_indentation < 0.0 {
            return Err(NormalThermalError::Invalid(
                "normal interval reaches separation; shorten or locate the contact event",
            ));
        }
        let motion_trials = Cell::new(0usize);
        let contact_trials = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            poll()?;
            if x.len() != 1 || !x[0].is_finite() || x[0] < 0.0 {
                return Err(NormalThermalError::Invalid(
                    "normal-motion trial leaves the nonnegative indentation branch",
                ));
            }
            if motion_trials.get() >= controls.maximum_motion_trials {
                return Err(NormalThermalError::TrialBudget);
            }
            motion_trials.set(motion_trials.get() + 1);
            let candidate = self.contact.step_candidate(
                thermal,
                target_s,
                x[0],
                controls.contact,
                energy_guess_j,
                &resolve,
                || (*cancel.borrow_mut())(),
            )?;
            contact_trials.set(
                contact_trials
                    .get()
                    .checked_add(candidate.report().contact_trials)
                    .ok_or(NormalThermalError::TrialBudget)?,
            );
            let force = candidate.report().work.mean_normal_force_n;
            if force < 0.0 {
                return Err(NormalThermalError::Invalid(
                    "normal-motion branch requires compressive force",
                ));
            }
            let residual = finite(
                finite(x[0] - free_indentation)? + finite(0.5 * dt * dt * inverse_mass * force)?,
            )?;
            Ok((candidate, residual))
        };
        let problem =
            FallibleFiniteDifference::new(1, controls.difference_step_m, |x: &[f64]| {
                evaluate(x).map(|(_, r)| vec![r])
            })
            .map_err(difference_error)?;
        let construction =
            NewtonKrylovState::new(&problem, vec![indentation_guess_m], controls.newton);
        if let Some(e) = problem.take_failure() {
            return Err(difference_error(e));
        }
        let mut solver = construction.map_err(NormalThermalError::Newton)?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
        if let Some(e) = problem.take_failure() {
            return Err(difference_error(e));
        }
        poll()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(NormalThermalError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(NormalThermalError::NotConverged(result.report));
        }
        let (candidate, _) = evaluate(&solver.x)?;
        let force = candidate.report().work.mean_normal_force_n;
        let impulse = finite(dt * force)?;
        let mut velocities = [0.0; 2];
        let mut displacements = [0.0; 2];
        let mut momentum = [0.0; 2];
        let mut translations = self.translations_m;
        let mut kinetic_change = 0.0;
        let mut external_work = 0.0;
        for body in 0..2 {
            let signed_impulse = if body == 0 { -impulse } else { impulse };
            let dv =
                finite(finite(finite(dt * body_forces_n[body])? + signed_impulse)? / masses[body])?;
            velocities[body] = finite(self.velocities_mps[body] + dv)?;
            let actual_dv = finite(velocities[body] - self.velocities_mps[body])?;
            let midpoint = finite(self.velocities_mps[body] + 0.5 * actual_dv)?;
            displacements[body] = finite(dt * midpoint)?;
            translations[body] = finite(translations[body] + displacements[body])?;
            momentum[body] = finite(
                finite(masses[body] * actual_dv)?
                    - finite(dt * body_forces_n[body])?
                    - signed_impulse,
            )?;
            kinetic_change = finite(kinetic_change + finite(masses[body] * actual_dv * midpoint)?)?;
            external_work =
                finite(external_work + finite(body_forces_n[body] * displacements[body])?)?;
        }
        let actual_delta = finite(displacements[0] - displacements[1])?;
        let indentation_residual =
            finite(candidate.contact().endpoint().approach_m - old_indentation - actual_delta)?;
        let contact_body_work = finite(-force * actual_delta)?;
        let report = candidate.report();
        let mechanical_defect =
            finite(kinetic_change + report.work.work_into_contact_j - external_work)?;
        let native_work_defect = finite(contact_body_work + report.work.work_into_contact_j)?;
        let total_defect = finite(
            kinetic_change
                + report.work.internal_energy_change_j
                + report.thermal.stored_energy_change_j
                - external_work
                - report.external_heat_j,
        )?;
        if momentum
            .iter()
            .any(|v| v.abs() > controls.momentum_tolerance_n_s)
            || indentation_residual.abs() > controls.displacement_tolerance_m
            || [mechanical_defect, native_work_defect, total_defect]
                .iter()
                .any(|v| v.abs() > controls.total_energy_tolerance_j)
        {
            return Err(NormalThermalError::Invalid(
                "normal body momentum, midpoint geometry or complete energy gate failed",
            ));
        }
        let (thermal, contact, contact_report) = candidate.into_parts();
        let inventory = self
            .inventory
            .each_ref()
            .map(|old| thermal.cells()[&old.id].clone());
        let pair = Self {
            contact,
            inventory,
            velocities_mps: velocities,
            translations_m: translations,
        };
        pair.kinetic_energy_j().map_err(model_error)?;
        poll()?;
        Ok(NormalMotionCandidate {
            thermal,
            pair,
            report: NormalMotionReport {
                contact: contact_report,
                nonlinear: result.report,
                motion_trials: motion_trials.get(),
                contact_trials: contact_trials.get(),
                body_displacements_m: displacements,
                relative_midpoint_velocity_mps: finite(actual_delta / dt)?,
                compressive_impulse_n_s: impulse,
                body_momentum_residual_n_s: momentum,
                indentation_residual_m: indentation_residual,
                kinetic_energy_change_j: kinetic_change,
                external_work_j: external_work,
                mechanical_energy_residual_j: mechanical_defect,
                total_energy_residual_j: total_defect,
            },
        })
    }
}
fn whole_bodies<E>(
    thermal: &DepositionSimulation<'_>,
    bodies: &[DepositedCell; 2],
) -> Result<(), NormalThermalError<E>> {
    if bodies[0].body_id == bodies[1].body_id
        || thermal.cells().values().any(|cell| {
            bodies
                .iter()
                .any(|body| cell.body_id == body.body_id && cell.id != body.id)
        })
    {
        return Err(NormalThermalError::Invalid(
            "each normally moving body must own exactly one complete caloric cell",
        ));
    }
    Ok(())
}
