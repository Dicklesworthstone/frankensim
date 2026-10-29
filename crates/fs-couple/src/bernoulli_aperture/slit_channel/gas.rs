//! Complete source-gas changes in the original no-slip momentum field.

use super::{AcousticRealizeError, FixedChannel, GasState, SlitChannel, invalid};

/// Native kinetic storage transfer at one fixed physical time.
/// The prescribing gas owner owes this signed work separately from any
/// background thermodynamic/mass/species balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChannelGasUpdate {
    /// Resolved kinetic storage before the intervention [J].
    pub kinetic_energy_before_j: f64,
    /// Resolved kinetic storage after the intervention [J].
    pub kinetic_energy_after_j: f64,
    /// Signed external work into this channel, equal to E_after - E_before [J].
    /// This is neither pressure-port work nor physical/numerical dissipation.
    pub parameter_work_j: f64,
}

impl ChannelGasUpdate {
    pub(super) fn unchanged(energy_j: f64) -> Self {
        Self {
            kinetic_energy_before_j: energy_j,
            kinetic_energy_after_j: energy_j,
            parameter_work_j: 0.0,
        }
    }
}

pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AcousticRealizeError> {
    if cancelled() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}

pub(super) fn admit(
    gas: &GasState,
    maximum_absolute_parameter_work_j: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AcousticRealizeError> {
    poll(cancelled)?;
    if !maximum_absolute_parameter_work_j.is_finite() || maximum_absolute_parameter_work_j < 0.0 {
        return Err(invalid(
            "channel gas transfer requires finite nonnegative parameter-work allowance",
        ));
    }
    gas.validate_coherence()
        .map_err(|_| invalid("channel gas transfer requires one coherent complete gas"))
}

impl SlitChannel {
    /// Replace the complete source gas at fixed geometry and retained time.
    ///
    /// Each transverse cell retains its axial momentum, so v_new equals
    /// rho_old/rho_new times v_old. Source viscosity supplies the new original
    /// no-slip matrix, and the actual kinetic storage difference is separately
    /// budgeted parameter work. A viscosity-only update changes future wall
    /// diffusion without changing instantaneous momentum or kinetic storage.
    ///
    /// No sample, pressure work or heat is fabricated. `last_step()` remains
    /// the original historical exchange; `stored_energy_j()` reads the new
    /// endpoint. This is a prescribed material intervention, not a mass,
    /// species or background heat solve. The original cell count, time step,
    /// geometry and exit model remain fixed. Existing bounded coefficient
    /// assembly is bracketed by cancellation checks; transfer polls each 64
    /// cells. No accepted source or field changes until all checks succeed.
    pub fn rebind_gas_controlled(
        &mut self,
        gas: &GasState,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ChannelGasUpdate, AcousticRealizeError> {
        admit(gas, maximum_absolute_parameter_work_j, cancelled)?;
        if gas == &self.gas {
            poll(cancelled)?;
            return Ok(ChannelGasUpdate::unchanged(self.stored_energy_j()));
        }
        let mut candidate = Self::from_gas(gas, self.spec)?;
        let update = candidate.profile.transfer_material_history(
            &self.profile,
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        poll(cancelled)?;
        *self = candidate;
        Ok(update)
    }
}

impl FixedChannel {
    /// Native momentum-preserving gas transfer on the actual section geometry.
    /// Uses the same atomic source/work semantics as `SlitChannel` and
    /// `CircularChannel`; no hydraulic-diameter or lumped-inertance replacement.
    pub fn rebind_gas_controlled(
        &mut self,
        gas: &GasState,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ChannelGasUpdate, AcousticRealizeError> {
        match self {
            Self::Slit(channel) => {
                channel.rebind_gas_controlled(gas, maximum_absolute_parameter_work_j, cancelled)
            }
            Self::Circular(channel) => {
                channel.rebind_gas_controlled(gas, maximum_absolute_parameter_work_j, cancelled)
            }
        }
    }
}
