//! Instantaneous generalized momentum transfer through the native frame mass.

use super::{Frame2d, Rayleigh, SolidError, mass::MassMatrix};
use fs_exec::reduce::Compensated;

const RELATIVE_ROUNDOFF: f64 = 512.0 * f64::EPSILON;

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}

/// One accepted instantaneous frame event. Entries use the original node-major
/// coordinates `(ux, uy, theta)` and the same fixed-reference inertia as Newmark.
/// These are generalized moments, not a claim of geometrically exact world
/// angular momentum for a frame undergoing finite rotation.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameImpulse {
    /// Resolved `v_after - v_before`: m/s, m/s and rad/s at each node.
    pub velocity_change: Vec<f64>,
    /// `M * velocity_change`: kg m/s, kg m/s and kg m²/s at each node.
    /// Supported entries include the inertial cross terms of consistent mass.
    pub generalized_momentum_change: Vec<f64>,
    /// Support impulse ON the frame, `M * velocity_change - applied_impulse`
    /// at fixed coordinates, zero at free coordinates. Fixed supports do no work.
    pub support_impulse: Vec<f64>,
    /// Native `1/2 v_before^T M v_before` [J].
    pub kinetic_energy_before_j: f64,
    /// Native `1/2 v_after^T M v_after` [J].
    pub kinetic_energy_after_j: f64,
    /// Stable kinetic change `v_mid^T M (v_after - v_before)` [J]. This avoids
    /// subtracting two large absolute energies when the event is small.
    pub kinetic_energy_change_j: f64,
    /// Signed work of the applied impulse, `v_mid^T applied_impulse` [J].
    /// Negative work removes kinetic energy; no strain/material work is added.
    pub impulse_work_j: f64,
}

/// Check the actual represented response, not just the algebraic solve. A
/// nonzero jump rounded away on a large stored velocity must not be accepted.
/// Each row has its own dimensioned allowance, with no absolute mass/force floor.
fn checked_action(
    mass: &MassMatrix,
    response: &[f64],
    target: &[f64],
    fixed: &[bool],
    what: &str,
) -> Result<Vec<f64>, SolidError> {
    if response.iter().any(|v| !v.is_finite()) {
        return Err(invalid(what));
    }
    let actual: Vec<f64> = (0..response.len())
        .map(|i| mass.row_action(i, response))
        .collect();
    for i in 0..response.len() {
        if !actual[i].is_finite() {
            return Err(invalid(what));
        }
        if !fixed[i] {
            let tolerance = RELATIVE_ROUNDOFF * target[i].abs()
                + mass.row_action_roundoff(i, response, RELATIVE_ROUNDOFF);
            if !tolerance.is_finite() || (actual[i] - target[i]).abs() > tolerance {
                return Err(invalid(what));
            }
        }
    }
    Ok(actual)
}

impl Frame2d {
    /// Apply one instantaneous generalized impulse, without advancing time.
    /// Translational entries have units N s; rotational entries have units N m s.
    /// The original free mass solves `M_ff * delta_v = impulse_f`; fixed
    /// coordinates retain zero velocity and receive the balancing support impulse.
    /// Displacements, section deformations and constitutive history do not change.
    ///
    /// `post_event_load` is the COMPLETE finite generalized right-hand side
    /// immediately after the event, in N and N m. A caller using base-relative
    /// motion includes its `-M r a_ground` term here. Together with the declared
    /// Rayleigh damping and original material viscosity this determines the
    /// right-limit finite acceleration `M a = load - C v_after - f_internal(u)`.
    /// The acceleration getter therefore never exposes a fabricated zero or a
    /// stale pre-event rate. Existing Rayleigh-reference semantics are retained.
    ///
    /// `continuation_dt` is the actual subsequent integration interval, used
    /// ONLY for existing source-domain/Nyquist admission of material viscosity.
    /// It is not an impulse duration: neither the velocity jump nor the force
    /// balance uses it, and no finite-duration force approximation is made.
    /// Subsequent Newmark calls reevaluate acceleration for their supplied loads.
    ///
    /// # Errors
    /// Requires positive inertia on every FREE coordinate. Free massless
    /// coordinates need a separate DAE velocity/acceleration projection and
    /// refuse even if their applied impulse is zero; use consistent member
    /// mass, explicit section rotary inertia or appropriate support constraints.
    /// Fixed massless coordinates and all-fixed support-only events are admitted.
    /// Wrong-sized/nonfinite input, invalid continuation context, failed native
    /// mass/damping/material admission, or an unrepresentable actual momentum,
    /// energy or acceleration response leaves the COMPLETE frame unchanged.
    pub fn apply_impulse(
        &mut self,
        impulse: &[f64],
        post_event_load: &[f64],
        damping: Rayleigh,
        continuation_dt: f64,
    ) -> Result<FrameImpulse, SolidError> {
        self.require_untimed_materials()?;
        let n = self.ndof();
        if impulse.len() != n
            || post_event_load.len() != n
            || impulse
                .iter()
                .chain(post_event_load)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "frame impulse and finite post-event load need one finite value per DOF",
            ));
        }
        if self.u.iter().chain(&self.vel).any(|v| !v.is_finite())
            || (0..n).any(|i| self.fixed[i] && self.vel[i] != 0.0)
        {
            return Err(invalid(
                "frame impulse requires finite accepted state and stationary fixed supports",
            ));
        }
        let mass = self.reference_mass()?;
        if (0..n).any(|i| !self.fixed[i] && mass.entry(i, i) <= 0.0) {
            return Err(invalid(
                "free massless coordinates require a DAE impulse projection; use consistent reference or explicit rotary inertia",
            ));
        }
        // Factor once for both the impulsive jump and its finite-force endpoint.
        let inverse = mass.prepare_solve(&self.fixed)?;
        let proposed = inverse.apply(impulse);
        let velocity: Vec<f64> = self
            .vel
            .iter()
            .zip(&proposed)
            .map(|(&old, &change)| if change == 0.0 { old } else { old + change })
            .collect();
        let velocity_change: Vec<f64> =
            velocity.iter().zip(&self.vel).map(|(a, b)| a - b).collect();
        let generalized_momentum_change = checked_action(
            &mass,
            &velocity_change,
            impulse,
            &self.fixed,
            "frame impulse cannot be represented by the accepted velocity change",
        )?;
        let support_impulse: Vec<f64> = (0..n)
            .map(|i| {
                if self.fixed[i] {
                    generalized_momentum_change[i] - impulse[i]
                } else {
                    0.0
                }
            })
            .collect();
        if support_impulse.iter().any(|v| !v.is_finite()) {
            return Err(invalid("frame support impulse is unrepresentable"));
        }

        let mut kinetic_change = Compensated::zero();
        let mut work = Compensated::zero();
        let mut work_allowance = Compensated::zero();
        for i in 0..n {
            let midpoint = 0.5 * self.vel[i] + 0.5 * velocity[i];
            let change = midpoint * generalized_momentum_change[i];
            let applied = midpoint * impulse[i];
            kinetic_change = kinetic_change.accumulate(change);
            work = work.accumulate(applied);
            work_allowance = work_allowance
                .accumulate(RELATIVE_ROUNDOFF * change.abs())
                .accumulate(RELATIVE_ROUNDOFF * applied.abs());
        }
        let kinetic_energy_before_j = mass.kinetic_energy(&self.vel);
        let kinetic_energy_after_j = mass.kinetic_energy(&velocity);
        let kinetic_energy_change_j = kinetic_change.value();
        let impulse_work_j = work.value();
        let work_allowance = work_allowance.value();
        if [
            kinetic_energy_before_j,
            kinetic_energy_after_j,
            kinetic_energy_change_j,
            impulse_work_j,
            work_allowance,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || kinetic_energy_before_j < 0.0
            || kinetic_energy_after_j < 0.0
            || (kinetic_energy_change_j - impulse_work_j).abs() > work_allowance
        {
            return Err(invalid(
                "frame impulse kinetic energy or work is unrepresentable",
            ));
        }

        let preparation = self.prepare_damping(&mass, damping, continuation_dt)?;
        let (internal_force, _, _) = self.assemble(&self.u)?;
        let finite_force: Vec<f64> = (0..n)
            .map(|i| {
                let damping_force: f64 = (0..n)
                    .map(|j| preparation.cmat[i * n + j] * velocity[j])
                    .sum();
                post_event_load[i] - damping_force - internal_force[i]
            })
            .collect();
        if finite_force.iter().any(|v| !v.is_finite()) {
            return Err(invalid("finite post-impulse force is unrepresentable"));
        }
        let acceleration = inverse.apply(&finite_force);
        checked_action(
            &mass,
            &acceleration,
            &finite_force,
            &self.fixed,
            "finite post-impulse acceleration does not resolve the original mass/force balance",
        )?;

        // No fallible operation follows publication, and no constitutive trial
        // is committed: this event changes rates, not positions or material state.
        let result = FrameImpulse {
            velocity_change,
            generalized_momentum_change,
            support_impulse,
            kinetic_energy_before_j,
            kinetic_energy_after_j,
            kinetic_energy_change_j,
            impulse_work_j,
        };
        self.vel = velocity;
        self.acc = acceleration;
        if let Some(reference) = preparation.new_reference {
            self.rayleigh_reference_stiffness = Some(reference);
        }
        Ok(result)
    }
}
