//! Ideal shaft brakes with simultaneous engagement and retained sticking.
//!
//! A locked rotor has p_i = J_i a_i.omega and constant relative phase. Its
//! axial inertia therefore joins the carrier tensor; only unlocked momenta
//! remain in the native gyrostat's stored vector. This is a reduced mechanical
//! mode, not repeated projection of an unconstrained drift. Engagement is a
//! perfectly inelastic internal event. No finite torque capacity, friction
//! coefficient, actuator dynamics or thermal destination is inferred.

use super::{
    HeldTorque, InertiaTensor3, LieStepError, Mat3, RotorBank, RotorImpulse, RotorState, RotorStep,
    So3, Vec3, components, dd, dot_dd, increment, invalid, product, scalar, subtract3, vector,
    work_residual,
};

const ROUNDING: f64 = 8192.0 * f64::EPSILON;

/// Prepared subset of ideal rotor/carrier shaft locks.
///
/// Preparation adds each selected axial inertia once and caches the ORIGINAL
/// native tensor representation. Engagement, impulses and drift are O(N),
/// allocation-free, and use no new equation solver. The canonical state stays
/// a [`RotorState`]: releasing a lock means passing the accepted momenta and
/// phases unchanged to the original bank or another prepared lock selection.
/// Release itself applies no impulse, resets no momentum and consumes no work.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotorLock<const N: usize> {
    bank: RotorBank<N>,
    locked: [bool; N],
    inertia: InertiaTensor3,
}

/// Simultaneous perfectly inelastic internal engagement at fixed attitude.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotorLockEngagement<const N: usize> {
    /// Original impulse owner, including actual carrier recoil and shaft work.
    pub impulse: RotorImpulse<N>,
    /// Nonnegative mechanical energy dissipated by all engaging locks, J.
    /// No per-lock partition or finite thermal destination is asserted.
    pub dissipated_energy_j: f64,
    /// Actual kinetic increment plus engagement loss; admitted at local scale.
    pub work_energy_residual_j: f64,
}

/// Simultaneous unlocked-shaft/external impulse while retaining ideal locks.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LockedRotorImpulse<const N: usize> {
    /// Actual impulses; selected shaft entries are the internal lock reactions.
    pub impulse: RotorImpulse<N>,
    /// Sum of signed actual work at the unlocked shaft ports, J.
    pub drive_work_j: f64,
    /// Measured ideal-reaction work, admitted as a rounding residual, not heat.
    pub reaction_work_j: f64,
    /// Kinetic increment minus drive/external work, excluding no physical loss.
    pub work_energy_residual_j: f64,
}

/// Driven finite motion with the selected ideal shaft locks retained.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LockedRotorStep<const N: usize> {
    /// Complete original rotor-bank trajectory. Locked axial-momentum changes
    /// include the reactions during both impulses and free precession.
    pub motion: RotorStep<N>,
    /// Signed work of the unlocked drives, excluding ideal lock reactions, J.
    pub drive_work_j: f64,
    /// Measured endpoint-kick reaction work, admitted as rounding, not heat.
    /// The free-mode reconstruction residual remains in the energy balance.
    pub reaction_work_j: f64,
    /// Kinetic change minus drive/external work and the native drift defect, J.
    /// This measured rounding residual does not reclassify a reaction as a source.
    pub work_energy_residual_j: f64,
}

/// Torque-free motion in the retained lock mode.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LockedRotorDrift<const N: usize> {
    pub state_after: RotorState<N>,
    /// Actual integrated internal reactions on selected absolute rotor momenta.
    /// Unlocked entries are exactly zero. Attitude transport is already included.
    pub axial_momentum_change_nms: [f64; N],
    /// Complete original carrier-plus-rotor kinetic increment, J.
    pub kinetic_energy_change_j: f64,
    /// Measured native reduced-gyrostat kinetic defect; never physical damping.
    pub drift_energy_defect_j: f64,
    /// Original kinetic increment minus native drift defect. Actual canonical
    /// reconstruction must resolve this at the local momentum-transport scale.
    pub work_energy_residual_j: f64,
    pub angular_momentum_residual_world: Vec3,
}

impl<const N: usize> RotorBank<N> {
    /// Prepare a simultaneous ideal lock selection on the current rotor bank.
    ///
    /// Engagement is explicit through [`RotorLock::engage`]; preparation does
    /// not change a state. All-false selection retains the original tensor and
    /// free dynamics. Fixed body axes, inertias and balanced rotors are required.
    pub fn lock_mode(&self, locked: [bool; N]) -> Result<RotorLock<N>, LieStepError> {
        let inertia = if locked.iter().any(|&value| value) {
            let mut matrix = self.inertia.matrix().m;
            for row in 0..3 {
                for col in row..3 {
                    let mut entry = dd(matrix[3 * row + col]);
                    for (i, rotor) in self.rotors.iter().enumerate() {
                        if locked[i] {
                            let axis = components(rotor.axis_body);
                            entry = entry
                                + dd(rotor.axial_inertia_kg_m2) * dd(axis[row]) * dd(axis[col]);
                        }
                    }
                    matrix[3 * row + col] = entry.to_f64();
                    matrix[3 * col + row] = entry.to_f64();
                }
            }
            InertiaTensor3::new(Mat3 { m: matrix })?
        } else {
            self.inertia
        };
        Ok(RotorLock {
            bank: *self,
            locked,
            inertia,
        })
    }
}

impl<const N: usize> RotorLock<N> {
    /// Selected shafts, in the original bank inventory order.
    pub const fn selected(&self) -> &[bool; N] {
        &self.locked
    }

    /// Reduced tensor of this mode, including only the selected axial inertias.
    pub const fn retained_inertia(&self) -> InertiaTensor3 {
        self.inertia
    }

    /// Engage all selected locks together, preserving total body momentum,
    /// attitude and every relative phase. Unselected rotor momenta are unchanged.
    ///
    /// The full shared mobility determines the reactions; nonorthogonal rotors
    /// are not stopped independently. Loss is the positive quadratic energy of
    /// the actual velocity change, and must close the original kinetic change.
    /// A state already within native lock resolution is retained exactly.
    pub fn engage(&self, state: &RotorState<N>) -> Result<RotorLockEngagement<N>, LieStepError> {
        self.bank.admit_state(state)?;
        let impulse = if self.admit_locked(state).is_ok() {
            self.bank.apply_impulse(state, [0.0; N], Vec3::default())?
        } else {
            self.solve_impulse(state, [0.0; N], Vec3::default())?
        };
        self.admit_locked(&impulse.state_after)?;
        let dc = self.bank.carrier_momentum(&impulse.state_after)?
            - self.bank.carrier_momentum(state)?;
        let mut loss = dot_dd(self.bank.inertia.velocity(dc)?, dc) * dd(0.5);
        for i in 0..N {
            let delta = impulse.axial_momentum_change_nms[i];
            loss = loss + dd(0.5 * delta / self.bank.rotors[i].axial_inertia_kg_m2) * dd(delta);
        }
        let loss = scalar(loss.to_f64(), "lock engagement loss")?;
        let residual = scalar(
            (dd(impulse.kinetic_energy_change_j) + dd(loss)).to_f64(),
            "lock engagement work",
        )?;
        if loss < 0.0 {
            return Err(invalid("nonnegative lock engagement loss"));
        }
        admit_work(residual.abs(), impulse.kinetic_energy_change_j.abs() + loss)?;
        Ok(RotorLockEngagement {
            impulse,
            dissipated_energy_j: loss,
            work_energy_residual_j: residual,
        })
    }

    /// Apply finite unlocked motor and external BODY impulses simultaneously.
    ///
    /// Selected shaft inputs must be zero: the corresponding reactions are
    /// solved, not prescribed. Input and output must satisfy the selected locks.
    /// Ideal reaction work is admitted independently of signed drive work; no
    /// engagement loss is added to a lock that was already held. Attitude and
    /// phases are unchanged. The original rounded-transfer admission still owns
    /// carrier recoil, and every error leaves the input untouched.
    pub fn apply_impulse(
        &self,
        state: &RotorState<N>,
        unlocked_impulse_nms: [f64; N],
        external_impulse_body: Vec3,
    ) -> Result<LockedRotorImpulse<N>, LieStepError> {
        self.admit_locked(state)?;
        // A held lock already admitted at native resolution is not projected
        // again when every physical port is zero. Tiny roundoff reactions would
        // otherwise demand an unresolved transfer from an unchanged state.
        let impulse = if unlocked_impulse_nms.iter().all(|&value| value == 0.0)
            && external_impulse_body == Vec3::default()
        {
            self.bank
                .apply_impulse(state, [0.0; N], external_impulse_body)?
        } else {
            self.solve_impulse(state, unlocked_impulse_nms, external_impulse_body)?
        };
        self.admit_locked(&impulse.state_after)?;
        let mut drive = dd(0.0);
        let mut reaction = dd(0.0);
        let mut absolute_reaction = 0.0;
        let mut scale = impulse.external_work_j.abs() + impulse.kinetic_energy_change_j.abs();
        for i in 0..N {
            let work = impulse.shaft_work_j[i];
            if self.locked[i] {
                reaction = reaction + dd(work);
                absolute_reaction += work.abs();
            } else {
                drive = drive + dd(work);
                scale += work.abs();
            }
        }
        let drive = scalar(drive.to_f64(), "unlocked drive work")?;
        let reaction = scalar(reaction.to_f64(), "ideal lock reaction work")?;
        let residual = scalar(
            (dd(impulse.kinetic_energy_change_j) - dd(drive) - dd(impulse.external_work_j))
                .to_f64(),
            "locked impulse work",
        )?;
        admit_work(absolute_reaction + residual.abs(), scale)?;
        Ok(LockedRotorImpulse {
            impulse,
            drive_work_j: drive,
            reaction_work_j: reaction,
            work_energy_residual_j: residual,
        })
    }

    /// Advance held unlocked-shaft drives and a body- or space-held external torque.
    ///
    /// Two original constrained half-impulses surround one retained-mode drift.
    /// The selected axial inertias stay in the cached carrier tensor throughout
    /// precession; their reaction impulses are solved, never supplied as drives.
    /// The input must already satisfy the locks, and selected shaft torques must
    /// be exactly zero. Engagement and its inelastic loss remain explicit.
    ///
    /// A spatial torque is transformed at both actual endpoint attitudes.
    /// Locked relative phases remain bit-for-bit fixed. Unlocked phases use the
    /// same kick-endpoint rate quadrature as shaft work, so a representable held
    /// drive obeys work = torque times its phase increment to rounding precision.
    /// The drift's native energy defect and canonical reconstruction roundoff
    /// remain numerical errors, never brake heat or unreported drive work.
    ///
    /// This is a symmetric second-order fixed-inertia approximation; resolve
    /// gyroscopic motion and phase by timestep refinement. Signed nonzero steps
    /// permit reversal checks. There is no finite brake capacity, state-dependent
    /// drive law, implicit engagement, allocation, or new equation solver.
    /// With no selected locks this delegates to the original bank step exactly.
    /// All native finite-value, reaction, work and phase refusals remain atomic.
    pub fn step(
        &self,
        state: &RotorState<N>,
        unlocked_torque_nm: [f64; N],
        external_torque: HeldTorque,
        duration_s: f64,
    ) -> Result<LockedRotorStep<N>, LieStepError> {
        if !self.locked.iter().any(|&locked| locked) {
            return self.step_receipt(self.bank.step(
                state,
                unlocked_torque_nm,
                external_torque,
                duration_s,
            )?);
        }
        scalar(duration_s, "retained-mode driven duration")?;
        let half = 0.5 * duration_s;
        if duration_s == 0.0 || half == 0.0 {
            return Err(invalid("nonzero resolved retained-mode driven duration"));
        }
        let mut shaft_impulse = [0.0; N];
        for i in 0..N {
            scalar(unlocked_torque_nm[i], "unlocked shaft torque")?;
            if self.locked[i] && unlocked_torque_nm[i] != 0.0 {
                return Err(invalid("prescribed torque on a locked shaft"));
            }
            shaft_impulse[i] = product(unlocked_torque_nm[i], half, "unlocked torque impulse")?;
        }
        let impulse_at = |rotation: So3| -> Result<Vec3, LieStepError> {
            let torque = match external_torque {
                HeldTorque::Body(value) => vector(value, "external body torque")?,
                HeldTorque::Space(value) => rotation.inverse().rotate(value)?,
            };
            Ok(Vec3::new(
                product(torque.x, half, "external torque impulse")?,
                product(torque.y, half, "external torque impulse")?,
                product(torque.z, half, "external torque impulse")?,
            ))
        };
        let first = self.apply_impulse(state, shaft_impulse, impulse_at(state.rotation)?)?;
        // The composing step owns phase quadrature. Do not round or admit an
        // unused drift-only phase before evaluating its actual drive-work rate.
        let drift = self.drift_momenta(&first.impulse.state_after, duration_s)?;
        let drifted = drift.state_after;
        let second = self.apply_impulse(&drifted, shaft_impulse, impulse_at(drifted.rotation)?)?;
        let mut phase = state.relative_angle_rad;
        for i in 0..N {
            if !self.locked[i] {
                let rate = 0.5 * first.impulse.mean_relative_rate_rad_s[i]
                    + 0.5 * second.impulse.mean_relative_rate_rad_s[i];
                let angle = product(duration_s, rate, "driven unlocked phase increment")?;
                phase[i] = increment(phase[i], angle, "driven retained-mode relative phase")?;
            }
        }
        let mut motion = self.bank.finish_step(
            state,
            &first.impulse,
            &drifted,
            &second.impulse,
            duration_s,
            Some(phase),
        )?;
        // Keep the reduced native defect distinct from reconstructing the
        // original rotor/carrier kinetic form, as in the standalone lock drift.
        motion.drift_energy_defect_j = drift.drift_energy_defect_j;
        motion.work_energy_residual_j = work_residual(
            motion.kinetic_energy_change_j,
            &motion.shaft_work_j,
            motion.external_work_j,
            motion.drift_energy_defect_j,
        )?;
        self.step_receipt(motion)
    }

    fn step_receipt(&self, motion: RotorStep<N>) -> Result<LockedRotorStep<N>, LieStepError> {
        let mut drive = dd(0.0);
        let mut reaction = dd(0.0);
        for i in 0..N {
            if self.locked[i] {
                reaction = reaction + dd(motion.shaft_work_j[i]);
            } else {
                drive = drive + dd(motion.shaft_work_j[i]);
            }
        }
        let drive_work_j = scalar(drive.to_f64(), "retained-mode drive work")?;
        let reaction_work_j = scalar(reaction.to_f64(), "retained-mode reaction work")?;
        let work_energy_residual_j = scalar(
            (dd(motion.kinetic_energy_change_j)
                - dd(drive_work_j)
                - dd(motion.external_work_j)
                - dd(motion.drift_energy_defect_j))
            .to_f64(),
            "retained-mode driven work",
        )?;
        Ok(LockedRotorStep {
            motion,
            drive_work_j,
            reaction_work_j,
            work_energy_residual_j,
        })
    }

    /// Advance torque-free constrained motion with the original gyrostat map.
    ///
    /// Selected axial inertia travels with the carrier; selected p therefore
    /// changes during precession while its relative phase stays bit-for-bit
    /// fixed. Unlocked p stays fixed and its phase uses endpoint trapezoidal
    /// relative rates. This is the native second-order map, with measured energy
    /// defect and a phase/refinement requirement, not an exact-energy claim.
    /// Signed nonzero duration permits free-motion reversal checks. Finite but
    /// unresolved canonical reactions or phase increments refuse atomically.
    pub fn drift(
        &self,
        state: &RotorState<N>,
        duration_s: f64,
    ) -> Result<LockedRotorDrift<N>, LieStepError> {
        let mut drift = self.drift_momenta(state, duration_s)?;
        let before_rates = self.bank.relative_rates(state)?;
        let after_rates = self.bank.relative_rates(&drift.state_after)?;
        for i in 0..N {
            if !self.locked[i] {
                let delta = product(
                    duration_s,
                    0.5 * before_rates[i] + 0.5 * after_rates[i],
                    "unlocked phase increment",
                )?;
                drift.state_after.relative_angle_rad[i] = increment(
                    state.relative_angle_rad[i],
                    delta,
                    "retained-mode relative phase",
                )?;
            }
        }
        Ok(drift)
    }

    fn drift_momenta(
        &self,
        state: &RotorState<N>,
        duration_s: f64,
    ) -> Result<LockedRotorDrift<N>, LieStepError> {
        self.admit_locked(state)?;
        let stored = self.unlocked_momentum(state)?;
        let mut after = *state;
        (after.rotation, after.total_momentum_body) = self.inertia.momentum_step(
            state.rotation,
            state.total_momentum_body,
            stored,
            duration_s,
        )?;
        let omega = self.inertia.velocity(after.total_momentum_body - stored)?;
        // Transport the already-admitted lock state instead of projecting its
        // absolute momentum again. In particular, coaxial drift has no body
        // momentum change and must not invent a rounding-sized reaction.
        let delta_omega = self
            .inertia
            .velocity(after.total_momentum_body - state.total_momentum_body)?;
        for i in 0..N {
            if self.locked[i] {
                let reaction = product(
                    self.bank.rotors[i].axial_inertia_kg_m2,
                    self.bank.rotors[i].axis_body.dot(delta_omega),
                    "retained locked axial momentum change",
                )?;
                after.axial_momentum_nms[i] = increment(
                    state.axial_momentum_nms[i],
                    reaction,
                    "retained locked axial momentum",
                )?;
            }
        }
        self.admit_locked(&after)?;
        let change =
            core::array::from_fn(|i| after.axial_momentum_nms[i] - state.axial_momentum_nms[i]);
        self.bank.admit_recoil(
            self.bank.carrier_momentum(state)?,
            self.bank.carrier_momentum(&after)?,
            after.total_momentum_body - state.total_momentum_body,
            &change,
        )?;
        let mean_omega = self
            .inertia
            .velocity(state.total_momentum_body - stored)?
            .scale(0.5)
            + omega.scale(0.5);
        let native_change = scalar(
            dot_dd(
                mean_omega,
                after.total_momentum_body - state.total_momentum_body,
            )
            .to_f64(),
            "native retained-mode energy change",
        )?;
        let kinetic_change = self.bank.energy_change(state, &after)?;
        let residual = scalar(
            (dd(kinetic_change) - dd(native_change)).to_f64(),
            "retained-mode energy reconstruction",
        )?;
        admit_work(residual.abs(), self.transport_work_scale(state, &after)?)?;
        let momentum_residual = subtract3(
            after.rotation.rotate(after.total_momentum_body)?,
            state.rotation.rotate(state.total_momentum_body)?,
            Vec3::default(),
        );
        vector(momentum_residual, "retained-mode spatial momentum")?;
        Ok(LockedRotorDrift {
            state_after: after,
            axial_momentum_change_nms: change,
            kinetic_energy_change_j: kinetic_change,
            drift_energy_defect_j: native_change,
            work_energy_residual_j: residual,
            angular_momentum_residual_world: momentum_residual,
        })
    }

    fn solve_impulse(
        &self,
        state: &RotorState<N>,
        mut impulse: [f64; N],
        external: Vec3,
    ) -> Result<RotorImpulse<N>, LieStepError> {
        let rates = self.bank.relative_rates(state)?;
        let actual_external =
            super::super::super::torque_kick(state.total_momentum_body, external, 1.0)?
                - state.total_momentum_body;
        let mut rhs = components(actual_external).map(dd);
        for i in 0..N {
            scalar(impulse[i], "unlocked shaft impulse")?;
            let rotor = self.bank.rotors[i];
            let transfer = if self.locked[i] {
                if impulse[i] != 0.0 {
                    return Err(invalid("prescribed impulse on a locked shaft"));
                }
                product(
                    rotor.axial_inertia_kg_m2,
                    rates[i],
                    "lock velocity correction",
                )?
            } else {
                let after = increment(state.axial_momentum_nms[i], impulse[i], "unlocked impulse")?;
                -(after - state.axial_momentum_nms[i])
            };
            for (value, axis) in rhs.iter_mut().zip(components(rotor.axis_body)) {
                *value = *value + dd(axis) * dd(transfer);
            }
        }
        let delta_omega =
            self.inertia
                .velocity(Vec3::new(rhs[0].to_f64(), rhs[1].to_f64(), rhs[2].to_f64()))?;
        for i in 0..N {
            if self.locked[i] {
                let rotor = self.bank.rotors[i];
                impulse[i] = scalar(
                    (dd(rotor.axial_inertia_kg_m2)
                        * (dot_dd(rotor.axis_body, delta_omega) - dd(rates[i])))
                    .to_f64(),
                    "simultaneous lock reaction",
                )?;
            }
        }
        self.bank.apply_impulse(state, impulse, external)
    }

    fn unlocked_momentum(&self, state: &RotorState<N>) -> Result<Vec3, LieStepError> {
        self.bank.rotor_momentum(&core::array::from_fn(|i| {
            if self.locked[i] {
                0.0
            } else {
                state.axial_momentum_nms[i]
            }
        }))
    }

    fn admit_locked(&self, state: &RotorState<N>) -> Result<(), LieStepError> {
        let actual = self.bank.angular_velocity(state)?;
        let expected = self
            .inertia
            .velocity(state.total_momentum_body - self.unlocked_momentum(state)?)?;
        let scale = actual.norm().max(expected.norm());
        if !scale.is_finite() || (actual - expected).norm() > ROUNDING * scale {
            return Err(invalid("unresolved retained-mode carrier velocity"));
        }
        for i in 0..N {
            if self.locked[i] {
                let absolute =
                    state.axial_momentum_nms[i] / self.bank.rotors[i].axial_inertia_kg_m2;
                let carried = self.bank.rotors[i].axis_body.dot(actual);
                let residual = scalar(absolute - carried, "locked relative shaft rate")?;
                let scale = absolute.abs().max(carried.abs());
                if residual.abs() > ROUNDING * scale {
                    return Err(invalid("state does not satisfy selected shaft locks"));
                }
            }
        }
        Ok(())
    }

    fn transport_work_scale(
        &self,
        before: &RotorState<N>,
        after: &RotorState<N>,
    ) -> Result<f64, LieStepError> {
        let c0 = self.bank.carrier_momentum(before)?;
        let c1 = self.bank.carrier_momentum(after)?;
        let average =
            self.bank.inertia.velocity(c0)?.scale(0.5) + self.bank.inertia.velocity(c1)?.scale(0.5);
        let mut scale = components(average)
            .into_iter()
            .zip(components(c1 - c0))
            .map(|(rate, impulse)| (rate * impulse).abs())
            .sum::<f64>();
        for i in 0..N {
            let average = (0.5 * before.axial_momentum_nms[i] + 0.5 * after.axial_momentum_nms[i])
                / self.bank.rotors[i].axial_inertia_kg_m2;
            scale += (average * (after.axial_momentum_nms[i] - before.axial_momentum_nms[i])).abs();
        }
        scalar(scale, "retained-mode local transport work scale")
    }
}

fn admit_work(residual: f64, scale: f64) -> Result<(), LieStepError> {
    scalar(residual, "locked-rotor work residual")?;
    scalar(scale, "locked-rotor work scale")?;
    if residual > ROUNDING * scale {
        return Err(invalid("unresolved locked-rotor work"));
    }
    Ok(())
}
