//! Motional induction in bulk conductors feeding the original winding circuit.
//!
//! Stationary stranded profiles retain the reciprocal terminal relation
//! V = delta(lambda)/dt + R I. Every current trial advances the actual moving
//! nonlinear field from accepted history, so its induced flux and screening
//! change the solved load current. No independently fitted speed voltage is
//! appended to the circuit. The field's conductor-frame current supplies both
//! Joule heat and Lorentz work, and private winding-port work cancels once.

use super::{
    FieldAccount, FieldTrial, NonlinearWindingCandidate, NonlinearWindingCircuit,
    NonlinearWindingControls, NonlinearWindingError, Result, poll,
};
use crate::planar::eddy::nonlinear::moving::{
    NonlinearMovingEddyControls, NonlinearMovingEddyStep,
};

impl<'m> NonlinearWindingCircuit<'m> {
    /// Solve the original native circuit against moving, source-resolved bulk
    /// conductors and reversible nonlinear or stationary fixed-recoil magnets.
    ///
    /// Velocity is a complete held P1 nodal field [m/s]. The original field owner
    /// checks divergence, connectivity, floating-region current, magnetic Peclet
    /// and true nonlinear residual limits. Mesh, source laws, conductivity and
    /// material occupancy stay fixed; this operation does not transport an
    /// interface or permanent magnetization. Every stranded winding cell and
    /// every recoil-magnet cell must have zero velocity. Otherwise a fixed
    /// winding profile would omit its own motion-induced terminal work.
    ///
    /// The endpoint-current, backward-Euler field scheme is first order. The
    /// returned account closes delta(storage) + heat + numerical loss + signed
    /// mechanical work = external source work. Negative mechanical work is
    /// energy supplied by the imposed-motion owner, not negative Joule heat.
    /// Rejected trials publish no field, circuit, linkage, current or clock.
    pub fn step_with_eulerian_velocity_candidate(
        &self,
        target_s: f64,
        velocity_m_s: &[[f64; 2]],
        external_inputs: &[f64],
        controls: NonlinearWindingControls<NonlinearMovingEddyControls>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NonlinearWindingCandidate<'m, NonlinearMovingEddyStep>> {
        poll(cancelled)?;
        let mesh = self.field.problem().mesh();
        if velocity_m_s.len() != mesh.vertices().len() {
            return Err(NonlinearWindingError::Invalid(
                "complete finite nodal winding/bulk velocity required",
            ));
        }
        for velocity in velocity_m_s {
            poll(cancelled)?;
            if velocity.iter().any(|v| !v.is_finite()) {
                return Err(NonlinearWindingError::Invalid(
                    "complete finite nodal winding/bulk velocity required",
                ));
            }
        }
        for source in &self.conductors {
            poll(cancelled)?;
            if mesh.faces()[source.triangle]
                .iter()
                .any(|&vertex| velocity_m_s[vertex as usize] != [0.0, 0.0])
            {
                return Err(NonlinearWindingError::Invalid(
                    "fixed stranded winding cells require zero material velocity",
                ));
            }
        }
        self.step_coupled(
            target_s,
            None,
            self.conductors.clone(),
            external_inputs,
            controls,
            |j, cancelled| {
                let mut check = || cancelled();
                let candidate = self.field.step_with_eulerian_velocity_candidate(
                    target_s,
                    velocity_m_s,
                    j,
                    controls.field,
                    &mut check,
                )?;
                let (state, step) = candidate.into_parts();
                let account = FieldAccount {
                    stored_before_j: step.stored_before_j,
                    stored_after_j: step.stored_after_j,
                    heat_j: step.heat_j,
                    supplied_work_j: step.supplied_work_j,
                    numerical_dissipation_j: step.numerical_dissipation_j,
                    parameter_work_j: 0.0,
                    mechanical_work_j: step.mechanical_work_j,
                };
                Ok(FieldTrial {
                    state,
                    step,
                    account,
                })
            },
            cancelled,
        )
    }
}
