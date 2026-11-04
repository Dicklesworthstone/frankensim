//! Instantaneous momentum jumps in the retained mass-normalized basis.

use super::{
    ModalAcousticFrame, ModalAcousticState, ModalAcousticTimeError, ModalAcousticTimeModel,
    ModalAcousticWorkspace, ModalStepData, stage_modal_frame,
};
use fs_math::dd::Dd;

impl ModalAcousticTimeModel {
    /// Apply already resolved generalized impulses `[N s / sqrt(kg)]` at a
    /// fixed configuration, returning an owned diagnostic frame.
    ///
    /// Mass normalization gives `v+ = v- + J` and `q+ = q-`. No physical time
    /// elapses: elastic forces and viscous damping supply no impulse, and
    /// dissipated heat is exactly zero. Signed input work is
    /// `sum J_k * (v-_k + v+_k)/2`, including braking and extraction of energy.
    /// Endpoint means, products and the signed sum retain double-double tails
    /// so opposite modal works do not erase smaller physical contributions.
    /// The nominal sample period is not used.
    ///
    /// This accepts an impulse from another physics owner; it does not solve
    /// restitution, contact, support reactions or omitted modal content. A
    /// finite-duration strike should use the force trajectory steppers. The
    /// pressure remains the existing narrow-band modal observation, not a
    /// broadband acoustic impact solution.
    ///
    /// # Errors
    /// Wrong impulse count, nonfinite input, a nonzero impulse entirely lost
    /// to velocity rounding, or the existing state, energy, pressure and
    /// work-balance limits. Refusals leave all modes unchanged.
    pub fn apply_impulse(
        &mut self,
        generalized_impulse_n_s_per_sqrt_kg: &[f64],
    ) -> Result<ModalAcousticFrame, ModalAcousticTimeError> {
        let mut workspace = ModalAcousticWorkspace::new(self);
        self.apply_impulse_into(generalized_impulse_n_s_per_sqrt_kg, &mut workspace)?;
        Ok(workspace.frame)
    }

    /// Allocation-free form of [`Self::apply_impulse`] using existing modal
    /// workspace storage. No transition coefficients or transcendental
    /// functions are evaluated, even on the first call.
    pub fn apply_impulse_into<'w>(
        &mut self,
        generalized_impulse_n_s_per_sqrt_kg: &[f64],
        workspace: &'w mut ModalAcousticWorkspace,
    ) -> Result<&'w ModalAcousticFrame, ModalAcousticTimeError> {
        if generalized_impulse_n_s_per_sqrt_kg.len() != self.modes.len() {
            return Err(ModalAcousticTimeError::ImpulseCountMismatch {
                expected: self.modes.len(),
                found: generalized_impulse_n_s_per_sqrt_kg.len(),
            });
        }
        self.apply_impulse_with(
            |index| Ok(generalized_impulse_n_s_per_sqrt_kg[index]),
            workspace,
        )
    }

    /// Apply a signed physical port impulse `[N s]` through reciprocal
    /// mode-shape weights `[1/sqrt(kg)]` and return an owned frame.
    ///
    /// The same shapes project `J_k = shape_k * J_port` and reconstruct the
    /// port velocity `v_port = sum shape_k * v_k`. Consequently the port work
    /// equals `J_port * (v_port- + v_port+)/2`, and its instantaneous inverse
    /// effective mass is `sum shape_k^2`. Free coordinates participate through
    /// exactly the same mass normalization as elastic modes.
    /// A nonzero projected impulse that underflows to zero refuses.
    ///
    /// Assumptions and refusal semantics match [`Self::apply_impulse`].
    pub fn apply_port_impulse(
        &mut self,
        port_shapes_per_sqrt_kg: &[f64],
        impulse_n_s: f64,
    ) -> Result<ModalAcousticFrame, ModalAcousticTimeError> {
        let mut workspace = ModalAcousticWorkspace::new(self);
        self.apply_port_impulse_into(port_shapes_per_sqrt_kg, impulse_n_s, &mut workspace)?;
        Ok(workspace.frame)
    }

    /// Allocation-free form of [`Self::apply_port_impulse`]. Projection and
    /// candidate admission share one bounded pass over the retained modes.
    pub fn apply_port_impulse_into<'w>(
        &mut self,
        port_shapes_per_sqrt_kg: &[f64],
        impulse_n_s: f64,
        workspace: &'w mut ModalAcousticWorkspace,
    ) -> Result<&'w ModalAcousticFrame, ModalAcousticTimeError> {
        if port_shapes_per_sqrt_kg.len() != self.modes.len() {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "impulse port shapes must match the admitted mode count",
            });
        }
        if !impulse_n_s.is_finite() {
            return Err(ModalAcousticTimeError::InvalidInput {
                what: "physical port impulse must be finite",
            });
        }
        self.apply_impulse_with(
            |index| {
                let shape = port_shapes_per_sqrt_kg[index];
                let projected = shape * impulse_n_s;
                if projected == 0.0 && shape != 0.0 && impulse_n_s != 0.0 {
                    return Err(ModalAcousticTimeError::InvalidInput {
                        what: "modal impulse projection loses nonzero momentum",
                    });
                }
                Ok(projected)
            },
            workspace,
        )
    }

    fn apply_impulse_with<'w>(
        &mut self,
        impulse_at: impl Fn(usize) -> Result<f64, ModalAcousticTimeError>,
        workspace: &'w mut ModalAcousticWorkspace,
    ) -> Result<&'w ModalAcousticFrame, ModalAcousticTimeError> {
        stage_modal_frame(
            &self.modes,
            &self.states,
            self.budget,
            workspace,
            "modal impulse work and kinetic energy change disagree",
            "aggregate impulse work and kinetic energy change disagree",
            |index, _mode, state| {
                let impulse = impulse_at(index)?;
                if !impulse.is_finite() {
                    return Err(ModalAcousticTimeError::InvalidInput {
                        what: "generalized modal impulses must be finite",
                    });
                }
                let velocity = state.velocity_m_sqrt_kg_per_s + impulse;
                if impulse != 0.0 && velocity == state.velocity_m_sqrt_kg_per_s {
                    return Err(ModalAcousticTimeError::InvalidInput {
                        what: "modal impulse is below representable velocity resolution",
                    });
                }
                let mean = Dd::from_f64(0.5 * state.velocity_m_sqrt_kg_per_s)
                    + Dd::from_f64(0.5 * velocity);
                Ok(ModalStepData {
                    next: ModalAcousticState {
                        velocity_m_sqrt_kg_per_s: velocity,
                        ..state
                    },
                    // Preserve mean/product tails until signed works have
                    // combined, without subtracting stored modal energies.
                    work: Dd::from_f64(impulse) * mean,
                    work_allowance: 0.0,
                    loss: 0.0,
                    loss_allowance: 0.0,
                })
            },
        )?;
        self.states.copy_from_slice(&workspace.candidate);
        Ok(&workspace.frame)
    }
}
