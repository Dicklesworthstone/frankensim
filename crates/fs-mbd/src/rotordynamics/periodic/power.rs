//! Power of the original held-load rotor cycle, including its stiffness switches.
//!
//! Native affine maps still own motion. The shared exact quadratic owner caches
//! integrals of `f.v`, `v.C.v` and `v.K_skew.q` on those same linear equations.
//! The actual endpoint momenta/configurations independently determine kinetic
//! and symmetric-stiffness storage changes. Each coefficient switch, including
//! the period seam, contributes `q.(K_next-K_old).q/2` at unchanged q/v.
//! Work can be negative: active bearings and circulatory stiffness are never
//! relabeled as heat. Gyroscopic work vanishes algebraically because the
//! original model admits exactly skew G; its full force remains in motion.
//!
//! Warm advances use O(segments*n²) cached contractions. This does not infer
//! continuous modulation between the authored switches, instantaneous peak
//! power, a heat destination, positivity of K/C, or stability of a forced orbit.

use super::{
    MassSolve, PreparedRotorCycle, RotorBearingModel, RotorError, RotorFloquetBudget,
    RotorPeriodicOrbit, RotorPeriodicSegment, RotorPeriodicTrajectory, finite, poll,
    prepare_forced_cycle, restore_entry,
};
use crate::rotordynamics::transient::RotorTransientState;
use fs_math::dd::Dd;
use fs_phs::linear_exact::{DenseExponentialBudget, quadratic_flow_integral_with_budget};

/// One original coefficient interval and its outgoing stiffness switch.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorSegmentWork {
    pub applied_work_j: f64,
    /// Signed `integral v.C.v dt`; negative values represent active supply.
    pub velocity_force_work_j: f64,
    /// Signed `integral v.K_skew.q dt`; excluded from elastic storage.
    pub nonconservative_stiffness_work_j: f64,
    pub kinetic_energy_change_j: f64,
    pub elastic_energy_change_j: f64,
    /// Work supplied by changing K at this endpoint, with no velocity reset.
    pub stiffness_switch_work_j: f64,
    /// Storage change + C/K_skew work - supplied force work, before switching.
    pub balance_residual_j: f64,
}

/// A complete cycle with separately observed source, material and parameter work.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorCyclePower {
    pub segments: Vec<RotorSegmentWork>,
    pub applied_work_j: f64,
    pub velocity_force_work_j: f64,
    pub nonconservative_stiffness_work_j: f64,
    pub stiffness_switch_work_j: f64,
    /// Endpoint storage in the FIRST segment's coefficients on both sides.
    pub kinetic_energy_change_j: f64,
    pub elastic_energy_change_j: f64,
    /// Delta storage + C/K_skew work - force work - stiffness-switch work.
    pub balance_residual_j: f64,
    /// Sum of absolute segment residuals; cancellation cannot hide defects.
    pub absolute_segment_residual_j: f64,
    pub average_applied_power_w: f64,
    pub average_velocity_force_power_w: f64,
    pub average_nonconservative_stiffness_power_w: f64,
    pub average_stiffness_switch_power_w: f64,
}

#[derive(Debug)]
pub enum RotorCyclePowerError {
    Dynamics(RotorError),
    /// Each segment, accumulated absolute defects and the full-cycle residual
    /// independently satisfy the same caller-authored absolute allowance.
    Energy { residual_j: f64, tolerance_j: f64 },
}
impl From<RotorError> for RotorCyclePowerError {
    fn from(value: RotorError) -> Self {
        Self::Dynamics(value)
    }
}
impl core::fmt::Display for RotorCyclePowerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "periodic rotor power refused: {self:?}")
    }
}
impl std::error::Error for RotorCyclePowerError {}

struct SegmentIntegrals<'a> {
    model: &'a RotorBearingModel,
    source_scale: f64,
    forms: [Vec<f64>; 3],
}

/// Native motion and exact work integrals, borrowing the original material
/// coefficients. No independent inertia, integrator or mutable state is added.
pub struct PreparedRotorPowerCycle<'a> {
    motion: PreparedRotorCycle,
    integrals: Vec<SegmentIntegrals<'a>>,
}

/// Prepare held-load motion and its physical work channels together.
///
/// The original whole-query budget includes all native Floquet work, three
/// Gramian preparations per segment, retained integral matrices and workspace.
/// All forms use the SAME mass/time coordinates and normalized source used by
/// native motion. Refusal/cancellation publishes no partially prepared cycle.
pub fn prepare_forced_power_cycle<'a>(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'a>],
    held_forces: &[&[f64]],
    relative_mass_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorPowerCycle<'a>, RotorError> {
    poll(&mut cancelled)?;
    let Some(first) = segments.first() else {
        return Err(RotorError::Invalid("nonempty power cycle required"));
    };
    let n = first.model.dofs();
    let e = 2 * n + 1;
    let count = segments.len();
    if count > super::MAX_ROTOR_PERIODIC_SEGMENTS {
        return Err(RotorError::PeriodicBudget("power cycle segment limit"));
    }
    let extra_bytes = (24 * e * e + count * (3 * e * e + 32)) * 8;
    let extra_work = (e as u64).pow(3) * 256 * 3 * count as u64;
    let native_budget = RotorFloquetBudget {
        max_dense_bytes: budget
            .max_dense_bytes
            .checked_sub(extra_bytes)
            .ok_or(RotorError::PeriodicBudget("power integral bytes"))?,
        max_multiply_terms: budget
            .max_multiply_terms
            .checked_sub(extra_work)
            .ok_or(RotorError::PeriodicBudget("power integral work"))?,
        ..budget
    };
    let motion = prepare_forced_cycle(
        spin_rad_s,
        segments,
        held_forces,
        relative_mass_residual_tolerance,
        native_budget,
        &mut cancelled,
    )?;
    let period = motion.report.period_seconds;
    let mass = MassSolve::new(first.model, &mut cancelled)?;
    let gramian_budget = DenseExponentialBudget {
        max_dimension: e,
        max_dense_bytes: 24 * e * e * 8,
        max_multiply_terms: (e as u64).pow(3) * 256,
    };
    let mut integrals = Vec::with_capacity(count);
    for (segment, force) in segments.iter().zip(held_forces) {
        poll(&mut cancelled)?;
        let (generator, _) = mass.segment_generator(
            segment,
            spin_rad_s,
            period,
            relative_mass_residual_tolerance,
            &mut cancelled,
        )?;
        let (source, _) = mass.force_column(
            force,
            segment.duration_seconds,
            period,
            relative_mass_residual_tolerance,
            &mut cancelled,
        )?;
        let maximum = source.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
        let source_scale = if maximum == 0.0 { 1.0 } else { maximum };
        let d = e - 1;
        let mut augmented = vec![0.0; e * e];
        for i in 0..d {
            augmented[i * e..i * e + d]
                .copy_from_slice(&generator[i * d..(i + 1) * d]);
            augmented[i * e + d] =
                represented(Dd::from_f64(source[i]) / Dd::from_f64(source_scale))?;
        }
        let mut forms: [Vec<f64>; 3] = core::array::from_fn(|_| vec![0.0; e * e]);
        let h = Dd::from_f64(segment.duration_seconds);
        let half = Dd::from_f64(0.5);
        for i in 0..n {
            poll(&mut cancelled)?;
            let vi = Dd::from_f64(motion.state_scale[n + i]);
            let supply = represented(
                half * h * Dd::from_f64(force[i]) / vi / Dd::from_f64(source_scale),
            )?;
            forms[0][(n + i) * e + d] = supply;
            forms[0][d * e + n + i] = supply;
            for j in 0..n {
                let vj = Dd::from_f64(motion.state_scale[n + j]);
                forms[1][(n + i) * e + n + j] = represented(
                    h * Dd::from_f64(segment.model.damping()[i * n + j]) / vi / vj,
                )?;
                let skew = half
                    * (Dd::from_f64(segment.model.stiffness()[i * n + j])
                        - Dd::from_f64(segment.model.stiffness()[j * n + i]));
                let entry = represented(
                    half * h * skew / vi / Dd::from_f64(motion.state_scale[j]),
                )?;
                forms[2][(n + i) * e + j] = entry;
                forms[2][j * e + n + i] = entry;
            }
        }
        for form in &mut forms {
            poll(&mut cancelled)?;
            // Identically absent channels have identically zero integrals;
            // symmetric stiffness does not need a separate cold exponential.
            if form.iter().all(|value| *value == 0.0) {
                continue;
            }
            let flow =
                quadratic_flow_integral_with_budget(&augmented, form, e, 1.0, gramian_budget)
                    .map_err(RotorError::PeriodicExponential)?;
            poll(&mut cancelled)?;
            *form = flow.integral;
        }
        integrals.push(SegmentIntegrals {
            model: segment.model,
            source_scale,
            forms,
        });
    }
    poll(&mut cancelled)?;
    Ok(PreparedRotorPowerCycle { motion, integrals })
}

impl PreparedRotorPowerCycle<'_> {
    pub fn motion_cycle(&self) -> &PreparedRotorCycle {
        &self.motion
    }

    /// Return native accepted motion and independently gated physical work.
    /// The borrowed initial state is unchanged on every refusal.
    pub fn advance(
        &self,
        initial: &RotorTransientState,
        relative_residual_tolerance: f64,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(RotorPeriodicTrajectory, RotorCyclePower), RotorCyclePowerError> {
        admit_energy_tolerance(energy_tolerance_j)?;
        let trajectory =
            self.motion.advance(initial, relative_residual_tolerance, &mut cancelled)?;
        let power = self.account(&trajectory, energy_tolerance_j, &mut cancelled)?;
        Ok((trajectory, power))
    }

    /// Compute the original unique periodic particular orbit and its actual
    /// load/modulator consumption. The original Floquet report retains any
    /// instability; the existence of this forced orbit is not attraction.
    pub fn periodic_orbit(
        &self,
        relative_residual_tolerance: f64,
        maximum_resolvent_gain: f64,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(RotorPeriodicOrbit, RotorCyclePower), RotorCyclePowerError> {
        admit_energy_tolerance(energy_tolerance_j)?;
        let orbit = self.motion.periodic_orbit(
            relative_residual_tolerance,
            maximum_resolvent_gain,
            &mut cancelled,
        )?;
        let power = self.account(&orbit.trajectory, energy_tolerance_j, &mut cancelled)?;
        Ok((orbit, power))
    }

    fn account(
        &self,
        trajectory: &RotorPeriodicTrajectory,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RotorCyclePower, RotorCyclePowerError> {
        let n = self.motion.report.dofs;
        let e = 2 * n + 1;
        let mut totals = [Dd::ZERO; 4];
        let mut absolute = Dd::ZERO;
        let mut segments = Vec::with_capacity(self.integrals.len());
        for (index, cached) in self.integrals.iter().enumerate() {
            poll(cancelled)?;
            let before = &trajectory.boundary_states[index];
            let after = &trajectory.boundary_states[index + 1];
            let mut z = [0.0; 2 * super::MAX_ROTOR_DOFS + 1];
            for (i, (&value, &scale)) in before
                .displacement
                .iter()
                .chain(&before.velocity)
                .zip(&self.motion.state_scale)
                .enumerate()
            {
                z[i] = restore_entry(value, 1.0, scale)?;
            }
            z[e - 1] = cached.source_scale;
            let mut work = [Dd::ZERO; 3];
            for (result, form) in work.iter_mut().zip(&cached.forms) {
                for i in 0..e {
                    poll(cancelled)?;
                    for j in 0..e {
                        *result = *result
                            + Dd::from_f64(z[i])
                                * Dd::from_f64(form[i * e + j])
                                * Dd::from_f64(z[j]);
                    }
                }
            }
            let kinetic =
                storage_change(cached.model.mass(), &before.velocity, &after.velocity);
            let elastic = storage_change(
                cached.model.stiffness(),
                &before.displacement,
                &after.displacement,
            );
            let next = self.integrals[(index + 1) % self.integrals.len()].model;
            let switch =
                switch_work(cached.model.stiffness(), next.stiffness(), &after.displacement);
            let residual = kinetic + elastic + work[1] + work[2] - work[0];
            let residual_j = rounded(residual)?;
            check_energy(residual_j, tolerance)?;
            absolute = absolute + Dd::from_f64(residual_j.abs());
            for i in 0..3 {
                totals[i] = totals[i] + work[i];
            }
            totals[3] = totals[3] + switch;
            segments.push(RotorSegmentWork {
                applied_work_j: rounded(work[0])?,
                velocity_force_work_j: rounded(work[1])?,
                nonconservative_stiffness_work_j: rounded(work[2])?,
                kinetic_energy_change_j: rounded(kinetic)?,
                elastic_energy_change_j: rounded(elastic)?,
                stiffness_switch_work_j: rounded(switch)?,
                balance_residual_j: residual_j,
            });
        }
        let first = &trajectory.boundary_states[0];
        let last = &trajectory.boundary_states[trajectory.boundary_states.len() - 1];
        let model = self.integrals[0].model;
        let kinetic = storage_change(model.mass(), &first.velocity, &last.velocity);
        let elastic = storage_change(model.stiffness(), &first.displacement, &last.displacement);
        let residual = rounded(kinetic + elastic + totals[1] + totals[2] - totals[0] - totals[3])?;
        let absolute_segment_residual_j = rounded(absolute)?;
        check_energy(residual, tolerance)?;
        check_energy(absolute_segment_residual_j, tolerance)?;
        let period = Dd::from_f64(self.motion.report.period_seconds);
        poll(cancelled)?;
        Ok(RotorCyclePower {
            segments,
            applied_work_j: rounded(totals[0])?,
            velocity_force_work_j: rounded(totals[1])?,
            nonconservative_stiffness_work_j: rounded(totals[2])?,
            stiffness_switch_work_j: rounded(totals[3])?,
            kinetic_energy_change_j: rounded(kinetic)?,
            elastic_energy_change_j: rounded(elastic)?,
            balance_residual_j: residual,
            absolute_segment_residual_j,
            average_applied_power_w: rounded(totals[0] / period)?,
            average_velocity_force_power_w: rounded(totals[1] / period)?,
            average_nonconservative_stiffness_power_w: rounded(totals[2] / period)?,
            average_stiffness_switch_power_w: rounded(totals[3] / period)?,
        })
    }
}

// The symmetric quadratic form is storage even when K itself is not symmetric.
// Polarization resolves small changes beside a large unchanged configuration.
fn storage_change(matrix: &[f64], before: &[f64], after: &[f64]) -> Dd {
    let n = before.len();
    let half = Dd::from_f64(0.5);
    let mut result = Dd::ZERO;
    for i in 0..n {
        let delta = Dd::from_f64(after[i]) - Dd::from_f64(before[i]);
        for j in 0..n {
            let symmetric =
                half * (Dd::from_f64(matrix[i * n + j]) + Dd::from_f64(matrix[j * n + i]));
            result = result
                + half * delta * symmetric * (Dd::from_f64(after[j]) + Dd::from_f64(before[j]));
        }
    }
    result
}
fn switch_work(before: &[f64], after: &[f64], q: &[f64]) -> Dd {
    let mut result = Dd::ZERO;
    for i in 0..q.len() {
        for j in 0..q.len() {
            result = result
                + Dd::from_f64(0.5)
                    * Dd::from_f64(q[i])
                    * (Dd::from_f64(after[i * q.len() + j])
                        - Dd::from_f64(before[i * q.len() + j]))
                    * Dd::from_f64(q[j]);
        }
    }
    result
}
fn rounded(value: Dd) -> Result<f64, RotorError> {
    finite(value.hi)?;
    finite(value.lo)?;
    finite(value.to_f64())
}
fn represented(value: Dd) -> Result<f64, RotorError> {
    let result = rounded(value)?;
    if result == 0.0 && value != Dd::ZERO {
        Err(RotorError::Invalid("power coordinate scaling underflow"))
    } else {
        Ok(result)
    }
}
fn admit_energy_tolerance(tolerance: f64) -> Result<(), RotorError> {
    if !tolerance.is_finite() || tolerance < 0.0 {
        Err(RotorError::Invalid("finite nonnegative cycle energy allowance required"))
    } else {
        Ok(())
    }
}
fn check_energy(residual_j: f64, tolerance_j: f64) -> Result<(), RotorCyclePowerError> {
    if residual_j.abs() > tolerance_j {
        Err(RotorCyclePowerError::Energy { residual_j, tolerance_j })
    } else {
        Ok(())
    }
}
