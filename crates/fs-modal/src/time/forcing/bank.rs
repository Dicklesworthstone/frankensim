//! Physical-coefficient derivatives bound to the admitted runtime mechanics.

use super::super::{ModalAcousticTimeError, ModalAcousticTimeModel};
use super::{
    LinearForceParameterPullback, LinearForceSensitivityFrame, LinearForceTransition,
    ModalAcousticState, PreparedLinearForceSensitivity, invalid,
};

impl ModalAcousticTimeModel {
    /// Cold-prepare exact stiffness and physical damping-rate sensitivities.
    ///
    /// This also prepares the original nominal-clock linear-force cache. Every
    /// derivative uses that very transition, including its authoritative rate
    /// for free/soft coordinates and its original ratio arithmetic for legacy
    /// banks. Diagnostic damping ratios are never reconstructed into a new law.
    ///
    /// Immutable preparation is shared by clones and repeated calls return
    /// immediately. Changing damping rates or ratios replaces the mechanics
    /// kernel, so the changed bank must prepare its own sensitivities before a
    /// query. A clone retaining the old mechanics keeps its valid old cache.
    /// Call this outside a real-time callback; the two query methods allocate
    /// nothing and evaluate no transcendental functions.
    ///
    /// # Errors
    /// An original ramp-preparation refusal or unrepresentable derivative
    /// coefficients. No partial sensitivity cache or physical state is published.
    pub fn prepare_linear_forcing_sensitivities(&self) -> Result<(), ModalAcousticTimeError> {
        if self.transition_kernel.linear_sensitivities.get().is_some() {
            return Ok(());
        }
        self.prepare_linear_forcing()?;
        let transitions = self.transition_kernel.linear.get().ok_or_else(|| {
            invalid("linear-force sensitivity preparation requires the original ramp cache")
        })?;
        let prepared = transitions
            .iter()
            .map(LinearForceTransition::prepare_sensitivities)
            .collect::<Result<Vec<_>, _>>()?
            .into_boxed_slice();
        // Another clone may have prepared this identical immutable kernel.
        let _ = self.transition_kernel.linear_sensitivities.set(prepared);
        Ok(())
    }

    /// Preview one original ramp step and its physical coefficient derivatives.
    ///
    /// The endpoint state matches [`Self::advance_exact_linear_step`]. Partials
    /// use independent mass-normalized stiffness k = omega^2 [1/s^2] and viscous
    /// rate c [1/s], with initial state, forces, clock and basis held fixed.
    /// A dimensional dashpot coefficient C on a fixed modal mass m therefore
    /// contributes (dJ/dc)/m. Changing mass or normalization also changes state,
    /// force and observation coordinates and requires those separate chain rules.
    ///
    /// Like the original preview, this pure query does not advance accepted
    /// state or apply the bank's state/energy budgets. Prepare with
    /// [`Self::prepare_linear_forcing_sensitivities`] before entering a callback.
    ///
    /// # Errors
    /// Missing preparation for the current mechanics, an invalid mode index,
    /// or nonfinite/unrepresentable scalar inputs, motion or sensitivities.
    pub fn advance_exact_linear_step_with_sensitivities(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        force_start: f64,
        force_end: f64,
    ) -> Result<LinearForceSensitivityFrame, ModalAcousticTimeError> {
        self.linear_sensitivity_transition(mode_index)?
            .advance(state, force_start, force_end)
    }

    /// Pull an endpoint objective through one admitted mode's exact ramp step.
    ///
    /// Returns initial-state, both endpoint-force, stiffness and damping-rate
    /// cotangents in the same physical parameterization as
    /// [`Self::advance_exact_linear_step_with_sensitivities`]. A transient
    /// reverse sweep adds observation cotangents at their actual endpoints,
    /// propagates both state cotangents, and accumulates shared coefficient and
    /// force-sample contributions. The runtime and its accepted states are read
    /// only. Material/basis/initial-state/force dependencies need their own chain
    /// rule; the bank does not infer a dimensional mass from a mode shape.
    ///
    /// # Errors
    /// Missing preparation for the current mechanics, an invalid mode index,
    /// or any refusal from [`PreparedLinearForceSensitivity::pullback`].
    pub fn pullback_exact_linear_step(
        &self,
        mode_index: usize,
        state: ModalAcousticState,
        force_start: f64,
        force_end: f64,
        endpoint_displacement_cotangent: f64,
        endpoint_velocity_cotangent: f64,
    ) -> Result<LinearForceParameterPullback, ModalAcousticTimeError> {
        self.linear_sensitivity_transition(mode_index)?.pullback(
            state,
            force_start,
            force_end,
            endpoint_displacement_cotangent,
            endpoint_velocity_cotangent,
        )
    }

    fn linear_sensitivity_transition(
        &self,
        mode_index: usize,
    ) -> Result<&PreparedLinearForceSensitivity, ModalAcousticTimeError> {
        self.transition_kernel
            .linear_sensitivities
            .get()
            .ok_or_else(|| {
                invalid(
                    "linear-force sensitivity queries require prepare_linear_forcing_sensitivities for the current mechanics",
                )
            })?
            .get(mode_index)
            .ok_or_else(|| invalid("linear-force sensitivity mode index is out of range"))
    }
}
