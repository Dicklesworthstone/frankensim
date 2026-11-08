//! Commit a physical basis change to the live oscillator bank.
//!
//! The existing transfer owner admits the original physical operators and
//! projects displacement and canonical momentum. The existing time owner
//! admits the replacement clock, damping, observation, state and energy.
//! Neither sorted mode numbers nor pressure amplitudes identify physical state.

use super::{ModalAcousticTimeError, ModalAcousticTimeModel, Sum, modal_energy};
use crate::{
    ModePair,
    transfer::{
        ModalMassTransfer, ModalRigidTransferControls, ModalTransferError,
        transfer_modal_state_with_rigid_modes,
    },
};
use fs_exec::Cx;
use fs_math::det;
use fs_sparse::Csr;

/// Original physical data for an instantaneous change on one reference space.
/// Source shapes bind the current runtime coordinates; target shapes bind the
/// replacement runtime coordinates. Both must use their declared mass norm.
/// The parameter owner supplies the target spectrum, damping and observation
/// model; no material law, force remapping or acoustic transfer is inferred.
#[derive(Clone, Copy, Debug)]
pub struct ModalBasisChange<'a> {
    /// Original SPD mass before the change.
    pub previous_mass: &'a Csr,
    /// Original SPD mass after the change, on the same DOFs and inertial frame.
    pub target_mass: &'a Csr,
    /// Original symmetric target stiffness, including its true rigid nullspace.
    pub target_stiffness: &'a Csr,
    /// Shapes/eigenvalues in the CURRENT runtime order, before the change.
    pub previous_modes: &'a [ModePair],
    /// Shapes/eigenvalues in the REPLACEMENT runtime order, after the change.
    pub target_modes: &'a [ModePair],
    /// Old equilibrium minus new equilibrium in the common coordinates [m].
    pub equilibrium_offset_m: &'a [f64],
}

/// Accepted update with independently checked physical and runtime storage.
#[derive(Clone, Debug)]
pub struct ModalRuntimeTransfer {
    /// Original physical projection, signed parameter work, discarded state,
    /// and canonical-momentum residual. Discarded storage is not viscous heat.
    pub physical: ModalMassTransfer,
    /// Storage of the actual time owner before the update [J].
    pub previous_runtime_energy_j: f64,
    /// Storage of the actual replacement time owner after the update [J].
    pub retained_runtime_energy_j: f64,
    /// H_new - H_old - W_parameter + E_discarded [J], independently contracted
    /// from the runtime's frequencies and states. Bounded by energy_tolerance_j.
    pub runtime_energy_residual_j: f64,
}

/// A refused update changes no state, coefficient, clock or prepared cache.
#[derive(Clone, Debug, PartialEq)]
pub enum ModalRuntimeTransferError {
    /// Source binding or clock mismatch; the original runtime is unchanged.
    Invalid(&'static str),
    /// Original physical transfer admission, projection or cancellation refusal.
    Transfer(ModalTransferError),
    /// Original replacement runtime state/energy admission refusal.
    Runtime(ModalAcousticTimeError),
}

impl core::fmt::Display for ModalRuntimeTransferError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MODAL-RUNTIME-TRANSFER: {self:?}")
    }
}
impl std::error::Error for ModalRuntimeTransferError {}
impl From<ModalTransferError> for ModalRuntimeTransferError {
    fn from(value: ModalTransferError) -> Self {
        Self::Transfer(value)
    }
}
impl From<ModalAcousticTimeError> for ModalRuntimeTransferError {
    fn from(value: ModalAcousticTimeError) -> Self {
        Self::Runtime(value)
    }
}

fn bind_spectrum(
    runtime: &ModalAcousticTimeModel,
    modes: &[ModePair],
) -> Result<(), ModalRuntimeTransferError> {
    if runtime.modes.len() != modes.len() {
        return Err(ModalRuntimeTransferError::Invalid(
            "physical basis must match the complete runtime mode count",
        ));
    }
    for (mode, physical) in runtime.modes.iter().zip(modes) {
        let omega = mode.angular_frequency_rad_s;
        let expected = det::sqrt(physical.lambda);
        // Permit only square-root rounding, not a changed physical spectrum.
        // In particular, a tiny nonzero frequency cannot stand in for a rigid
        // coordinate, or vice versa. The physical owner rechecks K*phi=lambda*M*phi.
        if !physical.lambda.is_finite()
            || physical.lambda < 0.0
            || (omega == 0.0) != (physical.lambda == 0.0)
            || (omega - expected).abs() > 8.0 * f64::EPSILON * omega.max(expected)
        {
            return Err(ModalRuntimeTransferError::Invalid(
                "runtime frequencies must match their original physical eigenvalues",
            ));
        }
    }
    Ok(())
}

fn energy(runtime: &ModalAcousticTimeModel) -> f64 {
    let mut sum = Sum::default();
    for (&mode, &state) in runtime.modes.iter().zip(&runtime.states) {
        sum.add(modal_energy(mode, state));
    }
    sum.value()
}

impl ModalAcousticTimeModel {
    /// Transfer live vibration and atomically install a physically rebuilt bank.
    ///
    /// `target` is an admitted, at-rest runtime on the SAME sample clock. Its
    /// frequencies, physical damping rates, pressure transfers, state budgets
    /// and cached exact transitions become authoritative only after the complete
    /// physical transfer and runtime admission pass. Its mode count may differ.
    /// The current state is read here, not supplied as a potentially stale copy.
    /// Existing cloned models remain independent and retain their old mechanics.
    ///
    /// The original transfer owner preserves physical displacement and canonical
    /// momentum before explicitly budgeted truncation, including mode crossings,
    /// degenerate rotations, changed inertia and complete rigid subspaces. It
    /// reports signed parameter work separately from discarded represented energy;
    /// this zero-duration update produces no force work or viscous heat. The
    /// actual runtime energies additionally close the same physical account.
    ///
    /// Build `target` and call this method outside an allocation-free audio
    /// callback: sparse source qualification and projection allocate. Target
    /// transition caches are reused after commit. External force/receiver maps
    /// and any constitutive or radiation memory must be transferred by their
    /// owners before the enclosing coupled transaction commits.
    ///
    /// This is a fixed-frame instantaneous Hamiltonian parameter update, not a
    /// continuously changing-coefficient integrator, remeshing, support impulse,
    /// mass-flux model or thermodynamic heat law. The caller owns physical event
    /// timing and the parameter work. Cancellation or any refusal leaves `self`
    /// unchanged; no partial new basis, state or kernel becomes visible.
    pub fn transfer_basis(
        &mut self,
        mut target: Self,
        change: ModalBasisChange<'_>,
        controls: ModalRigidTransferControls,
        cx: &Cx<'_>,
    ) -> Result<ModalRuntimeTransfer, ModalRuntimeTransferError> {
        cx.checkpoint().map_err(|_| ModalTransferError::Cancelled)?;
        if target.sample_period_s.to_bits() != self.sample_period_s.to_bits() {
            return Err(ModalRuntimeTransferError::Invalid(
                "basis transfer requires the same physical sample clock",
            ));
        }
        if target.states.iter().any(|state| {
            state.displacement_m_sqrt_kg != 0.0 || state.velocity_m_sqrt_kg_per_s != 0.0
        }) {
            return Err(ModalRuntimeTransferError::Invalid(
                "replacement bank must be at rest; another live history cannot be discarded",
            ));
        }
        bind_spectrum(self, change.previous_modes)?;
        bind_spectrum(&target, change.target_modes)?;
        let physical = transfer_modal_state_with_rigid_modes(
            change.previous_mass,
            change.target_mass,
            change.target_stiffness,
            change.previous_modes,
            change.target_modes,
            &self.states,
            change.equilibrium_offset_m,
            controls,
            cx,
        )?;
        target.restore_states(&physical.transfer.states)?;
        let previous_runtime_energy_j = energy(self);
        let retained_runtime_energy_j = energy(&target);
        let mut residual = Sum::default();
        residual.add(retained_runtime_energy_j);
        residual.add(-previous_runtime_energy_j);
        residual.add(-physical.transfer.parameter_work_j);
        residual.add(physical.transfer.kinetic_residual_j);
        residual.add(physical.transfer.elastic_residual_j);
        let runtime_energy_residual_j = residual.value();
        let limit = controls.transfer.transfer.energy_tolerance_j;
        if !runtime_energy_residual_j.is_finite() || runtime_energy_residual_j.abs() > limit {
            return Err(ModalTransferError::Limit {
                quantity: "runtime transfer energy residual",
                value: runtime_energy_residual_j.abs(),
                limit,
            }
            .into());
        }
        cx.checkpoint().map_err(|_| ModalTransferError::Cancelled)?;
        *self = target;
        Ok(ModalRuntimeTransfer {
            physical,
            previous_runtime_energy_j,
            retained_runtime_energy_j,
            runtime_energy_residual_j,
        })
    }
}
