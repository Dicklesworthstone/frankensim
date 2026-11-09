//! Held generalized loads and complete periodic orbits with original memories.

use super::{
    CycleForcing, MaxwellDampingBranch, PreparedRotorMaxwellCycle, RotorError, RotorFloquetBudget,
    RotorMaxwellPeriodicTrajectory, RotorMaxwellState, RotorPeriodicSegment, affine, finite, poll,
    prepare, restore_entry, state_values,
};
use fs_math::dd::Dd;

#[cfg(test)]
mod tests;

/// A unique numerical periodic particular orbit in original [q, q_dot, sigma].
///
/// Its Floquet spectrum belongs to the same full hereditary plan. Existence
/// of this forced solution does not imply attraction or attained steady state.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorMaxwellPeriodicOrbit {
    /// One complete physical period with the separately advanced endpoint.
    pub trajectory: RotorMaxwellPeriodicTrajectory,
    /// Original rowwise affine fixed-point residual on the actual initial state.
    /// Its scale uses (P-I)x and b, preserving small force and memory effects.
    pub relative_fixed_point_residual: f64,
    /// Endpoint/start mismatch in the common mass/time/material coordinates,
    /// relative to the complete-state peak over the returned boundaries.
    pub relative_closure_residual: f64,
    /// Native estimate of ||(I-P)^-1||_1 (1+||P||_1) in the common balanced
    /// coordinates. It is neither a norm bound nor a physical condition number.
    pub resolvent_gain_estimate: f64,
}

/// Prepare full-state affine cycles under held original generalized loads.
///
/// Each supplied row is the actual force/torque held over its corresponding
/// segment, in the model's original coordinate order and power-conjugate units.
/// A force change is not an impulse. C/K may change at the same boundary, while
/// displacement, velocity and every original branch stress remain continuous.
/// Mass, gyroscopic matrix, spin and branch definitions remain fixed.
///
/// K is RELAXED; the branch stiffness must not also be included in K. The load
/// enters only the original mechanical momentum rows. The full coupled native
/// exponential generates the induced material stresses; no zero-memory reset
/// or fitted viscous damping is used. Every original 2n+b state survives.
///
/// The shared normalized held-source exponential adds one nonphysical source
/// coordinate, for at most 129 coordinates. That source is excluded from the
/// physical period map and Floquet spectrum. Admission includes all retained
/// affine maps and one complete boundary trajectory/orbit solve. Cancellation
/// is polled around each bounded exponential and throughout assembly/solves.
///
/// Smooth load waveforms require independent segment refinement. Motion and
/// independently reconstructed residuals are returned; no actuator work, heat
/// ledger, nonlinear bearing law, eigenvalue enclosure or attraction is inferred.
pub fn prepare_maxwell_forced_cycle(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'_>],
    branches: &[MaxwellDampingBranch],
    held_forces: &[&[f64]],
    relative_generator_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorMaxwellCycle, RotorError> {
    prepare(
        spin_rad_s,
        segments,
        branches,
        None,
        CycleForcing::Held(held_forces),
        relative_generator_residual_tolerance,
        budget,
        true,
        cancelled,
    )
}

impl PreparedRotorMaxwellCycle {
    /// Original [q, q_dot, sigma] offset b in x_next = P x + b.
    pub fn forcing_offset(&self) -> &[f64] {
        &self.flow.physical_offset
    }

    /// Fixed physical boundary times [s], from zero through one full period.
    pub fn boundary_times_seconds(&self) -> &[f64] {
        &self.flow.boundary_times
    }

    /// Solve the full original-memory affine periodic equation (I-P)x=b.
    ///
    /// The same native partial-pivot LU and resolvent admission used by the
    /// mechanical periodic owner act on ALL 2n+b balanced states. No material
    /// rate, dark branch, repeated relaxation time or mechanical direction is
    /// removed. Singular maps refuse even when forcing admits a nonunique
    /// periodic family; no phase is pinned or loss added to regularize them.
    ///
    /// The caller's finite gain limit must be >=1. A near-identity period map
    /// cannot appear resolved merely because cond(I-P) is small. The native
    /// gain is an estimate, not an enclosure. The independently advanced full
    /// physical endpoint must close, and the original fixed-point and cycle
    /// increment equations must each pass the supplied tolerance in (0,1).
    /// Unstable cycles retain their positive Floquet growth unchanged.
    pub fn periodic_orbit(
        &self,
        relative_residual_tolerance: f64,
        maximum_resolvent_gain: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMaxwellPeriodicOrbit, RotorError> {
        poll(&mut cancelled)?;
        if !relative_residual_tolerance.is_finite()
            || relative_residual_tolerance <= 0.0
            || relative_residual_tolerance >= 1.0
        {
            return Err(RotorError::Invalid(
                "periodic memory orbit tolerance must be in (0,1)",
            ));
        }
        let d = self.flow.state_scale.len();
        let (balanced, resolvent_gain_estimate) = affine::solve_periodic_initial(
            &self.flow.monodromy,
            &self.flow.offset,
            maximum_resolvent_gain,
            &mut cancelled,
        )?;
        let mut original = Vec::with_capacity(d);
        for (&value, &scale) in balanced.iter().zip(&self.flow.state_scale) {
            poll(&mut cancelled)?;
            original.push(restore_entry(finite(value)?, scale, 1.0)?);
        }
        let n = self.report.dofs;
        let initial = RotorMaxwellState {
            displacement: original[..n].to_vec(),
            velocity: original[n..2 * n].to_vec(),
            branch_stress: original[2 * n..].to_vec(),
        };
        let relative_fixed_point_residual = affine::affine_residual(
            &self.report.monodromy,
            &self.flow.physical_offset,
            &original,
            &original,
            |_| Ok(0.0),
            relative_residual_tolerance,
            &mut cancelled,
        )?;
        let trajectory = self.advance(&initial, relative_residual_tolerance, &mut cancelled)?;
        // Exactly static or nodal components may acquire O(eps) roundoff.
        // Closure uses the common full-state norm; independent rowwise
        // affine checks above still use the original force/increment scale.
        let mut peak = 0.0_f64;
        for state in &trajectory.boundary_states {
            poll(&mut cancelled)?;
            for value in self.flow.encode(
                &state_values(state, self.report.dofs, self.report.memory_states)?,
                &mut cancelled,
            )? {
                peak = peak.max(value.abs());
            }
        }
        let last = state_values(
            &trajectory.boundary_states[trajectory.boundary_states.len() - 1],
            self.report.dofs,
            self.report.memory_states,
        )?;
        let mut closure = 0.0_f64;
        for i in 0..d {
            poll(&mut cancelled)?;
            let difference = (Dd::from_f64(last[i]) - Dd::from_f64(original[i]))
                * Dd::from_f64(self.flow.state_scale[i]);
            closure = closure.max(finite(difference.to_f64().abs())?);
        }
        let relative_closure_residual = if peak == 0.0 {
            closure
        } else {
            finite(closure / peak)?
        };
        if relative_closure_residual > relative_residual_tolerance {
            return Err(RotorError::Residual {
                relative: relative_closure_residual,
                limit: relative_residual_tolerance,
            });
        }
        poll(&mut cancelled)?;
        Ok(RotorMaxwellPeriodicOrbit {
            trajectory,
            relative_fixed_point_residual,
            relative_closure_residual,
            resolvent_gain_estimate,
        })
    }
}
