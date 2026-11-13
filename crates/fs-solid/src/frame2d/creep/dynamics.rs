//! Single-interval backward-Euler motion of the original timed fiber members.
//!
//! The existing force-based trial, consistent tangent, reference mass and
//! backtracking own mechanics. Acceleration is the nonlinear unknown, so no
//! initial constitutive evaluation consumes an extra creep interval. All trials
//! and the independent material audit start from the same accepted histories.

use super::super::backtrack;
use super::{
    Frame2d, FrameCreepEnergy, FrameCreepError, FrameCreepTolerances, FrameFiberCreepEnergy,
    Geometry, SolidError, energy_gate, poll, positive, require,
};
use fs_exec::reduce::Compensated;
use std::cell::Cell;

/// Complete mechanical work of one first-order inertial creep step.
#[derive(Clone, Copy, Debug, Default)]
pub struct FrameCreepDynamicEnergy {
    /// Native fiber storage/heat and actual external/support work. Its two
    /// residuals subtract `inertial_work_j` before auditing material work.
    /// `numerical_dissipation_j` here contains only native material loss.
    pub material: FrameCreepEnergy,
    pub initial_kinetic_j: f64,
    pub final_kinetic_j: f64,
    /// Stable `average(v) dot M delta(v)`, without subtracting large inventories.
    pub kinetic_change_j: f64,
    /// Endpoint `M a dot delta(u)`, using the resolved kinematic increment.
    pub inertial_work_j: f64,
    /// `0.5 delta(v)^T M delta(v)`: time-integration loss, never material heat.
    pub inertial_numerical_dissipation_j: f64,
    /// Inertial work minus kinetic change and numerical inertia loss.
    pub inertial_work_residual_j: f64,
    /// External/support work minus every material and kinetic storage/loss term.
    pub balance_residual_j: f64,
}

impl FrameCreepDynamicEnergy {
    /// Physical creep/plastic heat plus signed reciprocal thermal exchange.
    /// Neither material nor inertial numerical dissipation enters this transfer.
    pub fn caloric_heat_j(self) -> f64 {
        self.material.caloric_heat_j()
    }
}

/// Immutable proposal; accepting `frame` advances motion and history together.
#[derive(Clone, Debug)]
pub struct FrameCreepDynamicStep {
    pub frame: Frame2d,
    pub duration_s: f64,
    /// Global Newton corrections; element/section trials retain their own bounds.
    pub iterations: usize,
    pub fibers: Vec<FrameFiberCreepEnergy>,
    pub energy: FrameCreepDynamicEnergy,
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}

impl Frame2d {
    /// Advance one absolute load endpoint with the original frame inertia.
    ///
    /// Every timed fiber must carry exactly `duration_s`; retime the native laws
    /// explicitly with [`Self::with_creep_laws`] before a different interval.
    /// The old acceleration only seeds Newton and does not age the material.
    /// Newton ends when every free residual is within its force or moment
    /// allowance and the work those residuals do over the step is within
    /// `energy_j`.
    /// Stationary supports, small linear geometry and the native Euler–Bernoulli
    /// fiber members are admitted. Attached viscosity, transverse shear and a
    /// retained Rayleigh reference require their own dynamic work owner.
    ///
    /// This is a first-order, dissipative time integrator. Refine the physical
    /// interval to resolve motion and bound its reported numerical inertia loss.
    /// Physical creep/plastic loss and signed reciprocal heat retain the native
    /// material partition. Failure or cancellation leaves `self` unchanged.
    /// A nonzero position or rate increment that rounds away is refused;
    /// partially rounded positions must still pass the native material-work gate.
    #[allow(clippy::too_many_lines)]
    pub fn advance_creep_backward_euler(
        &self,
        load: &[f64],
        duration_s: f64,
        tolerances: FrameCreepTolerances,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FrameCreepDynamicStep, FrameCreepError> {
        poll(&mut cancelled)?;
        let n = self.ndof();
        let position_weight = duration_s * duration_s;
        require(
            self.geometry == Geometry::Linear
                && !self.elements.is_empty()
                && positive(duration_s)
                && positive(position_weight)
                && self.u.len() == n
                && self.vel.len() == n
                && self.acc.len() == n
                && self.vel.iter().chain(&self.acc).all(|v| v.is_finite())
                && (0..n).all(|d| !self.fixed[d] || (self.vel[d] == 0.0 && self.acc[d] == 0.0)),
            "creep motion requires finite state, stationary supports, linear geometry and a representable physical interval",
        )?;
        require(
            self.rayleigh_reference_stiffness.is_none()
                && self
                    .elements
                    .iter()
                    .all(|e| e.transverse_shear.is_none() && e.viscous_section.is_none()),
            "creep motion does not own Rayleigh, section viscosity or transverse-shear work",
        )?;
        let (count, material_duration) = self.admit_creep_step(load, tolerances, &mut cancelled)?;
        require(
            duration_s == material_duration
                && positive(duration_s * tolerances.force_n)
                && positive(duration_s * tolerances.moment_n_m),
            "inertia and all creep fibers need the same duration and representable impulse allowances",
        )?;
        let mass = self.reference_mass().map_err(FrameCreepError::Mechanical)?;
        let allowance = |d: usize| {
            if d % 3 == 2 {
                tolerances.moment_n_m
            } else {
                tolerances.force_n
            }
        };
        let velocity = |a: &[f64]| -> Vec<f64> {
            (0..n)
                .map(|d| duration_s.mul_add(a[d], self.vel[d]))
                .collect()
        };
        let position = |a: &[f64]| -> Vec<f64> {
            (0..n)
                .map(|d| self.u[d] + duration_s * duration_s.mul_add(a[d], self.vel[d]))
                .collect()
        };
        let residual = |a: &[f64], force: &[f64]| -> Vec<f64> {
            (0..n)
                .map(|d| load[d] - mass.row_action(d, a) - force[d])
                .collect()
        };
        let norm = |r: &[f64]| -> f64 {
            (0..n)
                .filter(|&d| !self.fixed[d])
                .fold(0.0_f64, |value, d| value.max(r[d].abs() / allowance(d)))
        };
        let (iterations, acceleration, endpoint, endpoint_velocity, forces, trials) = {
            let cancelled_seen = Cell::new(false);
            let mut checkpoint = || {
                if cancelled_seen.get() || cancelled() {
                    cancelled_seen.set(true);
                    Err(invalid("frame creep motion cancelled"))
                } else {
                    Ok(())
                }
            };
            let solve = (|| {
                let assemble = |a: &[f64], checkpoint: &mut _| {
                    let u = position(a);
                    let v = velocity(a);
                    if a.iter().chain(&u).chain(&v).any(|x| !x.is_finite()) {
                        return Err(invalid("creep endpoint kinematics are unrepresentable"));
                    }
                    let asm = self.assemble_controlled(&u, checkpoint)?;
                    if asm.0.iter().chain(&asm.1).any(|x| !x.is_finite())
                        || residual(a, &asm.0).iter().any(|x| !x.is_finite())
                    {
                        return Err(invalid(
                            "creep endpoint force or tangent is unrepresentable",
                        ));
                    }
                    Ok(asm)
                };
                let mut acceleration = self.acc.clone();
                let mut assembled = assemble(&acceleration, &mut checkpoint);
                if assembled.is_err() {
                    // A predictor can leave the admitted material domain while
                    // the implicit solution is valid. This seed targets the old
                    // position; it never commits its material response.
                    for (d, value) in acceleration.iter_mut().enumerate() {
                        *value = -self.vel[d] / duration_s;
                    }
                    assembled = assemble(&acceleration, &mut checkpoint);
                }
                let (mut force, mut tangent, mut trials) = assembled?;
                let mut history = Vec::new();
                for iteration in 0..60 {
                    checkpoint()?;
                    let r = residual(&acceleration, &force);
                    let rn = norm(&r);
                    history.push(rn);
                    // The force allowance alone does not end the solve: the
                    // work the free residual does over the step enters the
                    // equilibrium-work gate below, and at the allowance it can
                    // exceed the energy tolerance (1e-5 N over 2.6 mm is 2.6e-8
                    // J against 1e-8 J). Iterate until both are met.
                    let rate = velocity(&acceleration);
                    let residual_work: f64 = (0..n)
                        .filter(|&d| !self.fixed[d])
                        .map(|d| (r[d] * (duration_s * rate[d])).abs())
                        .sum();
                    if rn <= 1.0 && residual_work <= tolerances.energy_j {
                        return Ok((
                            iteration,
                            acceleration.clone(),
                            position(&acceleration),
                            velocity(&acceleration),
                            force,
                            trials,
                        ));
                    }
                    for entry in &mut tangent {
                        *entry *= position_weight;
                    }
                    mass.add_scaled(&mut tangent, 1.0);
                    if tangent.iter().any(|x| !x.is_finite()) {
                        return Err(invalid("creep effective tangent is unrepresentable"));
                    }
                    let correction =
                        self.solve_free(&tangent, &r)
                            .ok_or_else(|| SolidError::NewtonStalled {
                                history: history.clone(),
                            })?;
                    let (candidate, asm) = backtrack(
                        &acceleration,
                        &correction,
                        rn,
                        &mut checkpoint,
                        |a, checkpoint| {
                            let asm = assemble(a, checkpoint)?;
                            let rn = norm(&residual(a, &asm.0));
                            Ok((asm, rn))
                        },
                    )?;
                    acceleration = candidate;
                    (force, tangent, trials) = asm;
                }
                Err(SolidError::NewtonStalled { history })
            })();
            solve.map_err(|error| {
                if cancelled_seen.get() {
                    FrameCreepError::Cancelled
                } else {
                    FrameCreepError::Mechanical(error)
                }
            })?
        };
        poll(&mut cancelled)?;
        let delta_v: Vec<f64> = endpoint_velocity
            .iter()
            .zip(&self.vel)
            .map(|(v, old)| v - old)
            .collect();
        require(
            delta_v.iter().all(|x| x.is_finite()),
            "creep velocity increment is unrepresentable",
        )?;
        let mut applied_work = Compensated::zero();
        let mut support_work = Compensated::zero();
        let mut inertial_work = Compensated::zero();
        let mut kinetic_change = Compensated::zero();
        let mut inertia_loss = Compensated::zero();
        for d in 0..n {
            poll(&mut cancelled)?;
            let delta_u = duration_s * endpoint_velocity[d];
            let resolved_delta_v = duration_s * acceleration[d];
            require(
                delta_u.is_finite()
                    && resolved_delta_v.is_finite()
                    && (delta_u == 0.0 || endpoint[d] != self.u[d])
                    && (resolved_delta_v == 0.0 || delta_v[d] != 0.0),
                "creep position or rate increment cannot be represented in accepted coordinates",
            )?;
            let inertia = mass.row_action(d, &acceleration);
            let impulse = mass.row_action(d, &delta_v);
            let reaction = inertia + forces[d] - load[d];
            let impulse_residual = if self.fixed[d] {
                impulse - duration_s * inertia
            } else {
                impulse - duration_s * (load[d] - forces[d])
            };
            if !reaction.is_finite()
                || !impulse_residual.is_finite()
                || (!self.fixed[d] && reaction.abs() > allowance(d))
                || impulse_residual.abs() > duration_s * allowance(d)
            {
                return Err(FrameCreepError::Equilibrium {
                    dof: d,
                    residual: if !self.fixed[d] && reaction.abs() > allowance(d) {
                        reaction
                    } else {
                        impulse_residual / duration_s
                    },
                    tolerance: allowance(d),
                });
            }
            applied_work = applied_work.accumulate(load[d] * delta_u);
            if self.fixed[d] {
                support_work = support_work.accumulate(reaction * delta_u);
            }
            inertial_work = inertial_work.accumulate(inertia * delta_u);
            kinetic_change = kinetic_change
                .accumulate((0.5 * self.vel[d] + 0.5 * endpoint_velocity[d]) * impulse);
            inertia_loss = inertia_loss.accumulate(0.5 * delta_v[d] * impulse);
        }
        let mut energy = FrameCreepDynamicEnergy {
            initial_kinetic_j: mass.kinetic_energy(&self.vel),
            final_kinetic_j: mass.kinetic_energy(&endpoint_velocity),
            kinetic_change_j: kinetic_change.value(),
            inertial_work_j: inertial_work.value(),
            inertial_numerical_dissipation_j: inertia_loss.value(),
            ..FrameCreepDynamicEnergy::default()
        };
        require(
            [
                energy.initial_kinetic_j,
                energy.final_kinetic_j,
                energy.inertial_numerical_dissipation_j,
            ]
            .into_iter()
            .all(|value| value.is_finite() && value >= 0.0),
            "creep motion requires finite nonnegative kinetic storage and numerical loss",
        )?;
        energy.inertial_work_residual_j = energy.inertial_work_j
            - energy.kinetic_change_j
            - energy.inertial_numerical_dissipation_j;
        energy_gate(energy.inertial_work_residual_j, tolerances.energy_j)?;
        let mut next = self.clone();
        next.u = endpoint;
        next.vel = endpoint_velocity;
        next.acc = acceleration;
        next.commit_trials(trials);
        let (fibers, material_work) = self.audit_creep_materials(
            &next,
            count,
            tolerances,
            &mut energy.material,
            &mut cancelled,
        )?;
        energy.material.applied_work_j = applied_work.value();
        energy.material.support_work_j = support_work.value();
        let work = energy.material.applied_work_j + energy.material.support_work_j
            - energy.inertial_work_j;
        energy.material.equilibrium_work_residual_j = work - material_work;
        energy.material.balance_residual_j = work
            - (energy.material.final_stored_j - energy.material.initial_stored_j)
            - energy.material.thermoelastic_storage_change_j
            - energy.material.physical_dissipation_j
            - energy.material.reversible_heat_j
            - energy.material.numerical_dissipation_j;
        energy.balance_residual_j =
            energy.material.balance_residual_j + energy.inertial_work_residual_j;
        for residual in [
            energy.material.equilibrium_work_residual_j,
            energy.material.balance_residual_j,
            energy.balance_residual_j,
        ] {
            energy_gate(residual, tolerances.energy_j)?;
        }
        poll(&mut cancelled)?;
        Ok(FrameCreepDynamicStep {
            frame: next,
            duration_s,
            iterations,
            fibers,
            energy,
        })
    }
}
