//! Retained microslip, finite heat and two-body planar inertia in one candidate.
//!
//! The native contact return map receives the ACTUAL midpoint relative travel.
//! Its endpoint wrench supplies equal and opposite body impulses. Midpoint
//! kinematics close kinetic work exactly; native backward-Euler spring loss
//! remains numerical loss and is never deposited as heat. Both temperatures,
//! source friction and all retained strains are solved from accepted history.
//!
//! Each designated caloric cell represents one complete translating body.
//! The normal patch/load and tangent frame remain fixed, with no normal work,
//! rotation or moving thermal footprint. Elastic sticking retains finite shear
//! compliance; it does not impose zero relative velocity. This endpoint contact
//! composition is first order and requires timestep refinement, including at
//! stopping/reversal. It does not claim event-exact Coulomb arrest.

use super::*;
use fs_mbd_material::ResolvedBody;
use fs_thermal::enthalpy::DepositedCell;
use std::sync::Arc;

/// Momentum solve around the original finite caloric/contact solver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertialContactSolve {
    pub contact: CaloricContactSolve,
    /// Each outer trial calls the original caloric/contact owner once.
    pub maximum_motion_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Directional finite-difference perturbation of contact impulse [N s].
    pub difference_step_n_s: f64,
    /// Euclidean norm of EACH independently reconstructed body residual [N s].
    pub impulse_tolerance_n_s: f64,
    pub total_energy_tolerance_j: f64,
    /// Both unknowns and residuals are physical tangent impulse [N s].
    pub newton: NewtonKrylovConfig,
}

/// Both complete finite body masses and the original accepted loaded contact.
#[derive(Clone)]
pub struct InertialContactPair<'a> {
    contact: CaloricPartialSlip<'a>,
    inventory: [DepositedCell; 2],
    material_bodies: Option<Arc<[ResolvedBody; 2]>>,
    frame: TangentFrame,
    velocities_mps: [[f64; 2]; 2],
    translations_m: [[f64; 2]; 2],
}

pub struct InertialContactReport<S> {
    pub contact: CaloricContactReport<S>,
    pub nonlinear: NewtonReport,
    pub motion_trials: usize,
    /// Includes successful inner solves discarded by the outer iteration.
    pub contact_trials: usize,
    /// Actual A/B midpoint displacements in the retained tangent frame [m].
    pub body_displacements_m: [[f64; 2]; 2],
    pub relative_midpoint_velocity_mps: [f64; 2],
    pub contact_impulse_on_a_n_s: [f64; 2],
    pub body_momentum_residual_n_s: [[f64; 2]; 2],
    pub kinetic_energy_change_j: f64,
    pub external_work_j: f64,
    /// Delta K + native work INTO contact - exterior mechanical work.
    pub mechanical_energy_residual_j: f64,
    /// Delta(K + U_contact + H_caloric) + D_numerical - W_external - Q_external.
    pub total_energy_residual_j: f64,
}

pub struct InertialContactCandidate<'a, 't, S> {
    thermal: DepositionSimulation<'t>,
    pair: InertialContactPair<'a>,
    report: InertialContactReport<S>,
}

impl<'a, 't, S> InertialContactCandidate<'a, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn pair(&self) -> &InertialContactPair<'a> {
        &self.pair
    }
    pub const fn report(&self) -> &InertialContactReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        InertialContactPair<'a>,
        InertialContactReport<S>,
    ) {
        (self.thermal, self.pair, self.report)
    }
}

impl<'a> InertialContactPair<'a> {
    /// Use the admitted finite cells' authored masses. Every cell is the whole,
    /// uniformly represented body, not an arbitrary heated fraction of its mass.
    /// Translations start at zero at attachment; velocities and loaded native
    /// strain/kinetic histories are retained without relaxation or reset.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        contact: CaloricPartialSlip<'a>,
        frame: TangentFrame,
        velocities_mps: [[f64; 2]; 2],
    ) -> Result<Self, CaloricContactError> {
        if thermal.time_s() != contact.time_s()
            || temperatures(thermal, contact.curves, contact.cells)? != contact.temperature_k()
            || velocities_mps.iter().flatten().any(|v| !v.is_finite())
            || contact.checkpoint.state().region_elastic_twist_rad() != [0.0; 2]
            || contact
                .checkpoint
                .state()
                .reference_frame()
                .is_some_and(|old| old != frame)
        {
            return Err(CaloricContactError::Invalid(
                "synchronized complete bodies and the original translational contact frame required",
            ));
        }
        let a = thermal.cells().get(&contact.cells[0]).unwrap().clone();
        let b = thermal.cells().get(&contact.cells[1]).unwrap().clone();
        validate_whole_bodies(thermal, [&a, &b])?;
        for cell in [&a, &b] {
            if !cell.mass_kg.is_finite() || cell.mass_kg <= 0.0 {
                return Err(CaloricContactError::Invalid(
                    "positive finite body mass required",
                ));
            }
        }
        let pair = Self {
            contact,
            inventory: [a, b],
            material_bodies: None,
            frame,
            velocities_mps,
            translations_m: [[0.0; 2]; 2],
        };
        pair.kinetic_energy_j()?;
        Ok(pair)
    }

    /// Retain original source-resolved mass/geometry, with no fresh density
    /// evaluation during heating. This requires the contact's original pinned
    /// material friction binding to identify the ordered body material cards.
    /// Each complete mass and volume must match its caloric cell to roundoff.
    pub fn from_material_bodies(
        thermal: &DepositionSimulation<'_>,
        contact: CaloricPartialSlip<'a>,
        bodies: [ResolvedBody; 2],
        frame: TangentFrame,
        velocities_mps: [[f64; 2]; 2],
    ) -> Result<Self, CaloricContactError> {
        let mut pair = Self::new(thermal, contact, frame, velocities_mps)?;
        let binding =
            pair.contact
                .material_friction_binding()
                .ok_or(CaloricContactError::Invalid(
                    "source masses require the original material contact binding",
                ))?;
        for (index, body) in bodies.iter().enumerate() {
            let identity = binding.point.materials[index].content_hash();
            if body
                .regions()
                .iter()
                .any(|r| r.material().card_identity() != identity)
            {
                return Err(CaloricContactError::Invalid(
                    "source body and ordered contact material must agree",
                ));
            }
            let volume = body.regions().iter().try_fold(0.0, |sum, region| {
                finite::<Infallible>(sum + region.geometry().volume_m3())
            })?;
            let cell = &pair.inventory[index];
            if !same_inventory(body.inertia().mass(), cell.mass_kg)
                || !same_inventory(volume, cell.geometry.volume_m3)
            {
                return Err(CaloricContactError::Invalid(
                    "source body mass and volume differ from complete caloric inventory",
                ));
            }
        }
        pair.material_bodies = Some(Arc::new(bodies));
        pair.kinetic_energy_j()?;
        Ok(pair)
    }

    pub const fn contact(&self) -> &CaloricPartialSlip<'a> {
        &self.contact
    }
    pub const fn velocities_mps(&self) -> [[f64; 2]; 2] {
        self.velocities_mps
    }
    pub const fn translations_m(&self) -> [[f64; 2]; 2] {
        self.translations_m
    }
    pub const fn frame(&self) -> TangentFrame {
        self.frame
    }
    pub fn material_bodies(&self) -> Option<&[ResolvedBody; 2]> {
        self.material_bodies.as_deref()
    }
    pub fn masses_kg(&self) -> [f64; 2] {
        self.material_bodies.as_deref().map_or_else(
            || self.inventory.each_ref().map(|c| c.mass_kg),
            |b| b.each_ref().map(|body| body.inertia().mass()),
        )
    }
    pub fn kinetic_energy_j(&self) -> Result<f64, CaloricContactError> {
        let mut energy = 0.0;
        for (mass, velocity) in self.masses_kg().into_iter().zip(self.velocities_mps) {
            for v in velocity {
                energy = finite(energy + finite(0.5 * mass * v * v)?)?;
            }
        }
        Ok(energy)
    }

    /// Constant exterior tangent forces drive both real inertias. All motion
    /// trials restart the same accepted contact and finite caloric network.
    /// The callback supplies only OTHER heat sources/conductances, as in the
    /// original contact owner. No contact heat or mechanical work is added twice.
    ///
    /// `contact_impulse_guess_n_s` is the impulse ON A; `energy_guess_j` follows
    /// the original sorted cell order. Each budget is explicit: at most
    /// maximum_motion_trials complete inner solves, each with its own contact
    /// trial budget. A failed inner source query is propagated without fallback.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        body_forces_n: [[f64; 2]; 2],
        ownership: &GeneralizedWorkOwnership,
        controls: InertialContactSolve,
        contact_impulse_guess_n_s: [f64; 2],
        energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<InertialContactCandidate<'a, 't, S>, CaloricContactError<E>> {
        let cancel = RefCell::new(cancelled);
        let poll = || {
            if (*cancel.borrow_mut())() {
                Err(CaloricContactError::Cancelled)
            } else {
                Ok(())
            }
        };
        poll()?;
        controls
            .newton
            .validate()
            .map_err(CaloricContactError::Newton)?;
        let dt = finite(target_s - self.contact.time_s())?;
        if thermal.time_s() != self.contact.time_s()
            || dt <= 0.0
            || controls.maximum_motion_trials == 0
            || controls.maximum_newton_iterations == 0
            || body_forces_n.iter().flatten().any(|v| !v.is_finite())
            || contact_impulse_guess_n_s.iter().any(|v| !v.is_finite())
            || [
                controls.impulse_tolerance_n_s,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(CaloricContactError::Invalid(
                "invalid inertial contact state, forces or controls",
            ));
        }
        for old in &self.inventory {
            let now = thermal
                .cells()
                .get(&old.id)
                .ok_or(CaloricContactError::Invalid(
                    "complete moving caloric body disappeared",
                ))?;
            if now.body_id != old.body_id
                || now.segment_id != old.segment_id
                || now.material_id != old.material_id
                || now.mass_kg != old.mass_kg
                || now.specific_h_j_kg != old.specific_h_j_kg
                || now.geometry != old.geometry
            {
                return Err(CaloricContactError::Invalid(
                    "moving body's mass, thermal history, identity and reference geometry changed",
                ));
            }
        }
        validate_whole_bodies(thermal, self.inventory.each_ref())?;
        let masses = self.masses_kg();
        let motion_trials = Cell::new(0usize);
        let contact_trials = Cell::new(0usize);
        let evaluate = |impulse: &[f64]| {
            poll()?;
            if impulse.len() != 2 || motion_trials.get() >= controls.maximum_motion_trials {
                return Err(CaloricContactError::TrialBudget);
            }
            motion_trials.set(motion_trials.get() + 1);
            let mut velocities = [[0.0; 2]; 2];
            let mut displacements = [[0.0; 2]; 2];
            let mut relative_midpoint = [0.0; 2];
            for body in 0..2 {
                for axis in 0..2 {
                    let signed_impulse = if body == 0 {
                        impulse[axis]
                    } else {
                        -impulse[axis]
                    };
                    let velocity_change = finite(
                        finite(finite(dt * body_forces_n[body][axis])? + signed_impulse)?
                            / masses[body],
                    )?;
                    velocities[body][axis] =
                        finite(self.velocities_mps[body][axis] + velocity_change)?;
                    let midpoint = finite(self.velocities_mps[body][axis] + 0.5 * velocity_change)?;
                    displacements[body][axis] = finite(dt * midpoint)?;
                    relative_midpoint[axis] = finite(
                        relative_midpoint[axis] + if body == 0 { midpoint } else { -midpoint },
                    )?;
                }
            }
            let speed = finite(det::hypot(relative_midpoint[0], relative_midpoint[1]))?;
            // The native kinematic carrier factors velocity into a nonnegative
            // speed scale and dimensionless direction. No body rotation is
            // inferred from its historically named rolling_speed_mps field.
            let kinematics = PartialSlipKinematics {
                creepage: if speed == 0.0 {
                    [0.0; 2]
                } else {
                    [relative_midpoint[0] / speed, relative_midpoint[1] / speed]
                },
                rolling_speed_mps: speed,
                torsional_spin_rad_per_s: 0.0,
                contact_spin_rad_per_s: 0.0,
                dt_s: dt,
            };
            let candidate = self.contact.step_candidate(
                thermal,
                target_s,
                self.frame,
                kinematics,
                ownership,
                controls.contact,
                energy_guess_j,
                &resolve,
                || (*cancel.borrow_mut())(),
            )?;
            contact_trials.set(
                contact_trials
                    .get()
                    .checked_add(candidate.report().contact_trials)
                    .ok_or(CaloricContactError::TrialBudget)?,
            );
            let force = candidate.report().contact.tangent_force_n;
            let residual = [
                finite(impulse[0] - finite(dt * force[0])?)?,
                finite(impulse[1] - finite(dt * force[1])?)?,
            ];
            Ok((
                candidate,
                velocities,
                displacements,
                relative_midpoint,
                residual,
            ))
        };
        let problem =
            FallibleFiniteDifference::new(2, controls.difference_step_n_s, |x: &[f64]| {
                evaluate(x).map(|(_, _, _, _, r)| r.to_vec())
            })
            .map_err(difference_error)?;
        let construction = NewtonKrylovState::new(
            &problem,
            contact_impulse_guess_n_s.to_vec(),
            controls.newton,
        );
        if let Some(e) = problem.take_failure() {
            return Err(difference_error(e));
        }
        let mut solver = construction.map_err(CaloricContactError::Newton)?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
        if let Some(e) = problem.take_failure() {
            return Err(difference_error(e));
        }
        poll()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(CaloricContactError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(CaloricContactError::NotConverged(result.report));
        }
        let (candidate, velocities, displacements, relative_midpoint, _) = evaluate(&solver.x)?;
        let force = candidate.report().contact.tangent_force_n;
        let impulse = [finite(dt * force[0])?, finite(dt * force[1])?];
        let mut momentum_residual = [[0.0; 2]; 2];
        let mut kinetic_change = 0.0;
        let mut external_work = 0.0;
        let mut contact_body_work = 0.0;
        let mut translations = self.translations_m;
        for body in 0..2 {
            for axis in 0..2 {
                let dv = finite(velocities[body][axis] - self.velocities_mps[body][axis])?;
                let signed_impulse = if body == 0 {
                    impulse[axis]
                } else {
                    -impulse[axis]
                };
                momentum_residual[body][axis] = finite(
                    finite(masses[body] * dv)?
                        - finite(dt * body_forces_n[body][axis])?
                        - signed_impulse,
                )?;
                // Product of velocity change and midpoint avoids subtracting
                // large kinetic energies in a common translating frame.
                let midpoint = finite(self.velocities_mps[body][axis] + 0.5 * dv)?;
                kinetic_change = finite(kinetic_change + finite(masses[body] * dv * midpoint)?)?;
                external_work = finite(
                    external_work + finite(body_forces_n[body][axis] * displacements[body][axis])?,
                )?;
                contact_body_work = finite(
                    contact_body_work
                        + finite(
                            if body == 0 { force[axis] } else { -force[axis] }
                                * displacements[body][axis],
                        )?,
                )?;
                translations[body][axis] =
                    finite(translations[body][axis] + displacements[body][axis])?;
            }
            if det::hypot(momentum_residual[body][0], momentum_residual[body][1])
                > controls.impulse_tolerance_n_s
            {
                return Err(CaloricContactError::Invalid(
                    "independent body momentum gate failed",
                ));
            }
        }
        let r = candidate.report();
        let mechanical_defect = finite(
            kinetic_change + r.contact.generalized_work.work_into_interface_j - external_work,
        )?;
        let native_work_defect =
            finite(contact_body_work + r.contact.generalized_work.work_into_interface_j)?;
        let total_defect = finite(
            kinetic_change
                + r.contact_internal_energy_change_j
                + r.thermal.stored_energy_change_j
                + r.contact.generalized_work.numerical_dissipation_j
                - external_work
                - r.external_heat_j,
        )?;
        if [mechanical_defect, native_work_defect, total_defect]
            .iter()
            .any(|v| v.abs() > controls.total_energy_tolerance_j)
        {
            return Err(CaloricContactError::Invalid(
                "actual inertial work and complete finite contact energy do not close",
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
            material_bodies: self.material_bodies.clone(),
            frame: self.frame,
            velocities_mps: velocities,
            translations_m: translations,
        };
        poll()?;
        Ok(InertialContactCandidate {
            thermal,
            pair,
            report: InertialContactReport {
                contact: contact_report,
                nonlinear: result.report,
                motion_trials: motion_trials.get(),
                contact_trials: contact_trials.get(),
                body_displacements_m: displacements,
                relative_midpoint_velocity_mps: relative_midpoint,
                contact_impulse_on_a_n_s: impulse,
                body_momentum_residual_n_s: momentum_residual,
                kinetic_energy_change_j: kinetic_change,
                external_work_j: external_work,
                mechanical_energy_residual_j: mechanical_defect,
                total_energy_residual_j: total_defect,
            },
        })
    }
}

fn same_inventory(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
}

fn validate_whole_bodies<E>(
    thermal: &DepositionSimulation<'_>,
    bodies: [&DepositedCell; 2],
) -> Result<(), CaloricContactError<E>> {
    if bodies[0].body_id == bodies[1].body_id
        || thermal.cells().values().any(|cell| {
            bodies
                .iter()
                .any(|body| cell.body_id == body.body_id && cell.id != body.id)
        })
    {
        return Err(CaloricContactError::Invalid(
            "each moving body must own exactly one complete caloric cell",
        ));
    }
    Ok(())
}
