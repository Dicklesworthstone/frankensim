//! Instantaneous physical impulses through the original full native inertia.
//!
//! Point kinematics, constrained membership and mass are unchanged. This is an
//! event port for a caller-owned collision/actuator law, not a contact solver.

use super::{
    MovingLoadError, MovingShellMotion, ShellLoadPoint, ShellSpatialResultant, finite, invalid,
    linear, poll,
};

/// One integrated physical wrench ON the shell at its current accepted time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellPointImpulse {
    /// Reference material point, including its force lever arm.
    pub point: ShellLoadPoint,
    /// Cartesian linear impulse, N s.
    pub linear_impulse_ns: [f64; 3],
    /// Cartesian free angular impulse about that point, N m s.
    pub angular_impulse_nms: [f64; 3],
}

/// One accepted event, without elapsed time or a change in displacement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingShellImpulse {
    /// The unchanged accepted physical time, s.
    pub time_s: f64,
    /// Applied point impulses, including their reference-origin moments.
    pub applied_impulse: ShellSpatialResultant,
    /// Impulse ON the shell from its original pinned support constraints.
    pub support_impulse: ShellSpatialResultant,
    /// Momentum change from the actual velocity difference and full native mass,
    /// including attached payloads.
    pub momentum_change: ShellSpatialResultant,
    /// Actual change minus applied and support impulses.
    pub momentum_residual: ShellSpatialResultant,
    /// Applied impulse dotted with the actual mean material velocity, J.
    pub applied_work_j: f64,
    /// Support impulse dotted with the unchanged prescribed velocity, J.
    pub support_work_j: f64,
    /// `v_mean^T M delta_v`, using the actual accepted velocity difference.
    /// This avoids subtracting two large absolute kinetic energies.
    pub kinetic_energy_change_j: f64,
    /// Kinetic change minus applied and support work; strain is unchanged.
    pub energy_residual_j: f64,
}

impl MovingShellMotion<'_> {
    /// Apply localized event impulses to the original complete native inertia.
    ///
    /// The current mass factorization solves `M_ff delta_v_f = J_f`; fixed or
    /// driven support velocities stay continuous. Their reaction impulse is
    /// `M_cf delta_v_f - J_c`. Attached payload inertia, offsets and full
    /// rotational coupling remain in these original operators.
    ///
    /// Accepted velocity differences are recomputed AFTER floating addition.
    /// A free-row impulse lost beyond floating-point backward accuracy refuses
    /// rather than being reported as delivered. Impulse/work receipts use those
    /// actual differences and the held prescribed support velocities.
    ///
    /// Time, step count, displacement, support motion and endpoint forces stay
    /// unchanged. The original held-force equilibrium refreshes acceleration
    /// and regular support reactions before atomic publication. Repeated events
    /// reuse the initial mass LU; Newmark factors are unchanged. Viscous force
    /// supplies no impulse during this zero-duration event.
    ///
    /// The caller owns event timing, collision/actuator laws and avoiding
    /// duplicate application. No contact detection, restitution, constraint
    /// membership change or finite-rotation shell law is inferred.
    ///
    /// Invalid inputs, unrepresentable events or cancellation retain the full
    /// accepted state, held loads and support reactions.
    pub fn apply_impulses(
        &mut self,
        impulses: &[ShellPointImpulse],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MovingShellImpulse, MovingLoadError> {
        poll(&mut cancelled)?;
        let load = self.model.project_components(
            impulses.len(),
            impulses.iter().map(|impulse| {
                (
                    impulse.point,
                    impulse.linear_impulse_ns,
                    impulse.angular_impulse_nms,
                )
            }),
            &mut cancelled,
        )?;
        let mut increment = vec![0.0; self.state.velocity.len()];
        poll(&mut cancelled)?;
        self.mass_solver
            .solve_factored_into(&load.reduced, &mut increment)
            .map_err(linear)?;
        let mut next = self.state.clone();
        for (velocity, &change) in next.velocity.iter_mut().zip(&increment) {
            if change != 0.0 {
                *velocity += change;
            }
        }
        finite(&next.velocity)?;
        let delta: Vec<_> = next
            .velocity
            .iter()
            .zip(&self.state.velocity)
            .map(|(a, b)| a - b)
            .collect();
        finite(&delta)?;
        let mean: Vec<_> = self
            .state
            .velocity
            .iter()
            .zip(&next.velocity)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let full_delta = self.model.expand(&delta);
        let full_mean = self
            .model
            .expand_with_supports(&mean, self.supports.as_ref().map(|s| &s.velocity_m_s));
        let mass = &self.model.assembly.full_mass;
        let change = mass.apply(&full_delta);
        finite(&change)?;

        // Backward error is scaled by the actual row products, not a fixed
        // dimensional threshold or the much larger background velocity.
        let width = mass.dimension();
        for &row in &self.model.assembly.free_dofs {
            let scale = load.full[row].abs()
                + mass.values()[row * width..(row + 1) * width]
                    .iter()
                    .zip(&full_delta)
                    .map(|(m, v)| (m * v).abs())
                    .sum::<f64>();
            if !scale.is_finite()
                || (change[row] - load.full[row]).abs() > 4096.0 * f64::EPSILON * scale
            {
                return Err(invalid(
                    "shell event impulse is not representable in the accepted velocity",
                ));
            }
        }
        let support = |row: usize| {
            if self.model.free_index[row].is_none() {
                change[row] - load.full[row]
            } else {
                0.0
            }
        };
        let momentum_change = self.model.resultant(&change)?;
        let support_impulse = self.model.resultant_by(support)?;
        let momentum_residual = self
            .model
            .resultant_by(|i| change[i] - load.full[i] - support(i))?;
        let applied_work_j: f64 = load.full.iter().zip(&full_mean).map(|(j, v)| j * v).sum();
        let support_work_j: f64 = full_mean
            .iter()
            .enumerate()
            .map(|(i, v)| support(i) * v)
            .sum();
        let kinetic_energy_change_j: f64 = change.iter().zip(&full_mean).map(|(p, v)| p * v).sum();
        let energy_residual_j = kinetic_energy_change_j - applied_work_j - support_work_j;
        finite(&[
            applied_work_j,
            support_work_j,
            kinetic_energy_change_j,
            energy_residual_j,
        ])?;

        let held = self
            .model
            .assembly
            .free_dofs
            .iter()
            .map(|&i| self.full_load[i])
            .collect();
        let rhs = self.model.acceleration_rhs(
            &next.displacement,
            &next.velocity,
            held,
            self.supports.as_ref(),
        )?;
        poll(&mut cancelled)?;
        self.mass_solver
            .solve_factored_into(&rhs, &mut next.acceleration)
            .map_err(linear)?;
        finite(&next.acceleration)?;
        next.kinetic_energy_j = mass.quadratic_energy(&self.model.expand_with_supports(
            &next.velocity,
            self.supports.as_ref().map(|s| &s.velocity_m_s),
        ));
        next.total_energy_j = next.kinetic_energy_j + next.strain_energy_j;
        finite(&[next.kinetic_energy_j, next.total_energy_j])?;
        let reactions = self
            .model
            .reactions(&next, self.supports.as_ref(), &self.full_load)?;
        poll(&mut cancelled)?;
        self.state = next;
        self.reactions = reactions;
        Ok(MovingShellImpulse {
            time_s: self.state.time_s,
            applied_impulse: load.resultant,
            support_impulse,
            momentum_change,
            momentum_residual,
            applied_work_j,
            support_work_j,
            kinetic_energy_change_j,
            energy_residual_j,
        })
    }
}
