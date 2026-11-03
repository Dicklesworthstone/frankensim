//! Original-mass impulse balance with prescribed instantaneous support motion.

use super::*;
use crate::dynamics::LinearSolidImpulseReport;
use crate::dynamics::spectral::SpectralImpulseCandidate;
use fs_math::dd::Dd;

/// An instantaneous event at the existing accepted time. Full vectors use
/// `represented_dofs()` order; support reactions use `support_dofs()` order.
#[derive(Clone, Debug)]
pub struct SupportedSolidImpulse {
    pub time_s: f64,
    /// Actual original consistent-mass momentum increment [kg m/s].
    pub momentum_change_kg_m_s: Vec<f64>,
    /// Impulse exerted BY the prescribed supports ON the solid [kg m/s].
    pub support_reaction_kg_m_s: Vec<f64>,
    /// Applied plus support impulse minus represented M*Delta v [kg m/s].
    pub impulse_residual_kg_m_s: Vec<f64>,
    pub maximum_impulse_residual_kg_m_s: f64,
    /// Midpoint velocity dotted with original M*Delta v, avoiding subtraction
    /// of two large kinetic-energy totals [J].
    pub kinetic_change_j: f64,
    /// Inner complete-coordinate response to J_f - M_fc*Delta v_c. Its work
    /// belongs to the partitioned free problem; full physical work is below.
    pub free: LinearSolidImpulseReport,
    pub energy: SupportedSolidEnergy,
}

/// The full original-coordinate response remains private until surrounding
/// contact laws admit the actual represented endpoint. Prepared operators and
/// the accepted free trajectory are borrowed, not cloned.
pub(super) struct SupportedImpulseCandidate {
    pub(super) velocity_mps: Vec<f64>,
    pub(super) energy: SupportedSolidEnergy,
    free: SpectralImpulseCandidate,
    momentum_change_kg_m_s: Vec<f64>,
    support_reaction_kg_m_s: Vec<f64>,
    impulse_residual_kg_m_s: Vec<f64>,
    maximum_impulse_residual_kg_m_s: f64,
    kinetic_change_j: f64,
    applied_correction_j: f64,
    support_correction_j: f64,
}

fn represented(value: Dd) -> Result<f64> {
    let result = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && result.is_finite() {
        Ok(result)
    } else {
        Err(invalid("unrepresentable support impulse or work"))
    }
}

fn product_sum(a: &[f64], b: &[f64]) -> Result<f64> {
    represented(a.iter().zip(b).fold(Dd::ZERO, |sum, (&a, &b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    }))
}

impl SupportedSpectralSolid {
    /// Apply an impulse in ORIGINAL represented-coordinate order and prescribe
    /// the support velocities immediately AFTER the event [m/s]. Supports may
    /// retain their current velocities, or jump; positions and time stay fixed.
    /// Applied support-row impulses are allowed and enter actuator reactions.
    ///
    /// The original consistent mass gives J_eff = J_f - M_fc*Delta v_c and
    /// R_c = M_cf*Delta v_f + M_cc*Delta v_c - J_c. Neither cross mass nor
    /// support work is discarded. Signed impulse work uses average pre/post
    /// velocity, and both full and free energy histories retain this event.
    /// The explicit impulse allowance applies to every ORIGINAL row [kg m/s].
    ///
    /// One O(n) candidate borrows the prepared operators; no matrix clone or
    /// refactorization occurs. All original-mass, work and cancellation checks
    /// precede publication of either owner. Preparation, displacement, support
    /// acceleration and accepted-step count remain unchanged. The retained
    /// acceleration governs subsequent continuous intervals from the new rate.
    /// This supplies an impulse actuator, not a contact/restitution law.
    pub fn apply_momentum_impulse(
        &mut self,
        applied_impulse_kg_m_s: &[f64],
        support_velocity_after_mps: &[f64],
        maximum_impulse_error_kg_m_s: f64,
        cx: &Cx<'_>,
    ) -> Result<SupportedSolidImpulse> {
        let trial = self.trial_momentum_impulse(
            applied_impulse_kg_m_s,
            support_velocity_after_mps,
            maximum_impulse_error_kg_m_s,
            cx,
        )?;
        poll(cx)?;
        Ok(self.accept_momentum_impulse(trial))
    }

    pub(super) fn trial_momentum_impulse(
        &self,
        applied_impulse_kg_m_s: &[f64],
        support_velocity_after_mps: &[f64],
        maximum_impulse_error_kg_m_s: f64,
        cx: &Cx<'_>,
    ) -> Result<SupportedImpulseCandidate> {
        poll(cx)?;
        let n = self.coordinates.len();
        if applied_impulse_kg_m_s.len() != n
            || support_velocity_after_mps.len() != self.support_rows.len()
            || applied_impulse_kg_m_s
                .iter()
                .chain(support_velocity_after_mps)
                .any(|value| !value.is_finite())
            || !maximum_impulse_error_kg_m_s.is_finite()
            || maximum_impulse_error_kg_m_s < 0.0
        {
            return Err(invalid(
                "complete finite impulses, support velocities and nonnegative impulse allowance required",
            ));
        }
        let mut velocity = copy(&self.velocity)?;
        let mut delta = zeros(n)?;
        for (c, &row) in self.support_rows.iter().enumerate() {
            poll(cx)?;
            velocity[row] = support_velocity_after_mps[c];
            delta[row] = velocity[row] - self.velocity[row];
            if !delta[row].is_finite() {
                return Err(invalid("support velocity increment overflow"));
            }
        }
        let mut effective = zeros(self.free_rows.len())?;
        for (f, &row) in self.free_rows.iter().enumerate() {
            poll(cx)?;
            let mut impulse = Dd::from_f64(applied_impulse_kg_m_s[row]);
            for &column in &self.support_rows {
                impulse = impulse
                    - Dd::from_f64(self.mass[row * n + column]) * Dd::from_f64(delta[column]);
            }
            effective[f] = represented(impulse)?;
        }
        let trial =
            self.free
                .trial_momentum_impulse(&effective, maximum_impulse_error_kg_m_s, cx)?;
        for (f, &row) in self.free_rows.iter().enumerate() {
            velocity[row] = trial.velocity_mps[f];
            delta[row] = velocity[row] - self.velocity[row];
        }
        let mut momentum_change_kg_m_s = zeros(n)?;
        let mut mean = zeros(n)?;
        for row in 0..n {
            poll(cx)?;
            momentum_change_kg_m_s[row] = product_sum(&self.mass[row * n..(row + 1) * n], &delta)?;
            mean[row] = 0.5 * self.velocity[row] + 0.5 * velocity[row];
            if !mean[row].is_finite() {
                return Err(invalid("impulsive midpoint velocity overflow"));
            }
        }
        let mut support_reaction_kg_m_s = zeros(self.support_rows.len())?;
        let mut reaction_by_row = zeros(n)?;
        let mut support_work = Dd::ZERO;
        for (c, &row) in self.support_rows.iter().enumerate() {
            poll(cx)?;
            let reaction = represented(
                Dd::from_f64(momentum_change_kg_m_s[row])
                    - Dd::from_f64(applied_impulse_kg_m_s[row]),
            )?;
            support_reaction_kg_m_s[c] = reaction;
            reaction_by_row[row] = reaction;
            support_work = support_work + Dd::from_f64(reaction) * Dd::from_f64(mean[row]);
        }
        let mut impulse_residual_kg_m_s = zeros(n)?;
        let mut maximum_impulse_residual_kg_m_s = 0.0_f64;
        for row in 0..n {
            poll(cx)?;
            impulse_residual_kg_m_s[row] = represented(
                Dd::from_f64(applied_impulse_kg_m_s[row]) + Dd::from_f64(reaction_by_row[row])
                    - Dd::from_f64(momentum_change_kg_m_s[row]),
            )?;
            let defect = impulse_residual_kg_m_s[row].abs();
            if !defect.is_finite() || defect > maximum_impulse_error_kg_m_s {
                return Err(LinearSolidImpulseError::Momentum {
                    defect_kg_m_s: defect,
                    allowance_kg_m_s: maximum_impulse_error_kg_m_s,
                }
                .into());
            }
            maximum_impulse_residual_kg_m_s = maximum_impulse_residual_kg_m_s.max(defect);
        }
        let applied_work_j = product_sum(applied_impulse_kg_m_s, &mean)?;
        let support_work_j = represented(support_work)?;
        let kinetic_change_j = product_sum(&momentum_change_kg_m_s, &mean)?;
        let kinetic_after_j = quadratic(&self.mass, &velocity, cx)?;
        if kinetic_after_j < 0.0 || !(self.stored_j + kinetic_after_j).is_finite() {
            return Err(invalid("invalid post-impulse physical energy"));
        }
        let (mut cumulative_applied_work_j, mut applied_correction_j) =
            (self.applied_work_j, self.applied_correction_j);
        let (mut cumulative_support_work_j, mut support_correction_j) =
            (self.support_work_j, self.support_correction_j);
        if applied_work_j != 0.0 {
            add(
                &mut cumulative_applied_work_j,
                &mut applied_correction_j,
                applied_work_j,
            );
        }
        if support_work_j != 0.0 {
            add(
                &mut cumulative_support_work_j,
                &mut support_correction_j,
                support_work_j,
            );
        }
        if [
            cumulative_applied_work_j,
            applied_correction_j,
            cumulative_support_work_j,
            support_correction_j,
        ]
        .iter()
        .any(|value| !value.is_finite())
        {
            return Err(invalid("accumulated impulse work overflow"));
        }
        let energy = SupportedSolidEnergy {
            stored_before_j: self.stored_j,
            stored_after_j: self.stored_j,
            kinetic_before_j: self.kinetic_j,
            kinetic_after_j,
            applied_work_j,
            support_work_j,
            cumulative_applied_work_j,
            cumulative_support_work_j,
            balance_error_j: product_sum(&mean, &impulse_residual_kg_m_s)?,
            // Preparation and continuous motion define the energy reference
            // from the rounded sum of both stores. Retain that convention:
            // splitting this sum here would give a zero event a new datum.
            cumulative_balance_error_j: represented(
                Dd::from_f64(cumulative_applied_work_j) + Dd::from_f64(cumulative_support_work_j)
                    - Dd::from_f64(self.stored_j + kinetic_after_j)
                    + Dd::from_f64(self.initial_energy_j),
            )?,
        };
        for defect_j in [energy.balance_error_j, energy.cumulative_balance_error_j] {
            if defect_j.abs() > self.controls.maximum_work_error_j {
                return Err(SupportedSolidError::Energy {
                    defect_j,
                    allowance_j: self.controls.maximum_work_error_j,
                });
            }
        }
        poll(cx)?;
        Ok(SupportedImpulseCandidate {
            velocity_mps: velocity,
            energy,
            free: trial,
            momentum_change_kg_m_s,
            support_reaction_kg_m_s,
            impulse_residual_kg_m_s,
            maximum_impulse_residual_kg_m_s,
            kinetic_change_j,
            applied_correction_j,
            support_correction_j,
        })
    }

    pub(super) fn accept_momentum_impulse(
        &mut self,
        trial: SupportedImpulseCandidate,
    ) -> SupportedSolidImpulse {
        let free = self.free.accept_momentum_impulse(trial.free);
        self.velocity = trial.velocity_mps;
        self.kinetic_j = trial.energy.kinetic_after_j;
        self.applied_work_j = trial.energy.cumulative_applied_work_j;
        self.applied_correction_j = trial.applied_correction_j;
        self.support_work_j = trial.energy.cumulative_support_work_j;
        self.support_correction_j = trial.support_correction_j;
        SupportedSolidImpulse {
            time_s: free.time_s,
            momentum_change_kg_m_s: trial.momentum_change_kg_m_s,
            support_reaction_kg_m_s: trial.support_reaction_kg_m_s,
            impulse_residual_kg_m_s: trial.impulse_residual_kg_m_s,
            maximum_impulse_residual_kg_m_s: trial.maximum_impulse_residual_kg_m_s,
            kinetic_change_j: trial.kinetic_change_j,
            free,
            energy: trial.energy,
        }
    }
}
