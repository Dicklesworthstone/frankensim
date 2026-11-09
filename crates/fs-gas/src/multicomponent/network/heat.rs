//! Original implicit wall heat on prepared, unpublished compact inventories.
use super::Plenum;
use crate::multicomponent::wall_exchange::{
    WallContact, WallControl, WallDrag, WallExchangeError, WallThermal, WallUpdate, exchange_cell,
};
use crate::multicomponent::{MixtureGas, poll};
use fs_qty::Time;

impl Plenum {
    /// Prepare a stationary rigid-tank heat candidate in existing storage.
    ///
    /// The original gas/wall source owns its implicit caloric solve, physical
    /// heat and finite-wall energy account. Prescribed-temperature reservoirs
    /// and constant-capacity finite walls use their existing source types.
    /// A caloric wall without its original material-law binding refuses.
    ///
    /// Every canonical species density, the volume and zero bulk momentum are
    /// retained. This adds no reaction heat, port enthalpy or pressure work.
    /// The returned contact contains the candidate finite-wall temperature;
    /// the composing owner must publish it together with the output tank.
    ///
    /// The output needs the exact source species count. All source checks and
    /// the last cancellation poll precede output publication. Warm calls do
    /// not allocate. This operation changes no clock. Original backward-Euler
    /// heat integration is first order even within symmetric source splitting.
    #[allow(clippy::too_many_arguments)]
    pub fn heat_candidate_into(
        &self,
        gas: &MixtureGas,
        duration: Time,
        thermal: WallThermal,
        control: WallControl,
        output: &mut Self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WallUpdate, WallExchangeError> {
        poll(&mut cancelled)?;
        if output.state.partial_densities.len() != gas.species_count() {
            return Err(WallExchangeError::Invalid(
                "complete mixture plenum heat output species axis required",
            ));
        }
        let update = exchange_cell(
            gas,
            self.state(),
            self.volume_m3,
            duration,
            WallContact {
                drag: WallDrag::None,
                thermal,
            },
            control,
            &mut cancelled,
        )?;
        if update.momentum != 0.0
            || update.report.wall_impulse_ns != 0.0
            || update.report.mechanical_work_to_gas_j != 0.0
            || !(self.volume_m3 * update.energy).is_finite()
        {
            return Err(WallExchangeError::Invalid(
                "unrepresentable stationary plenum heat inventory",
            ));
        }
        poll(&mut cancelled)?;
        output
            .state
            .partial_densities
            .copy_from_slice(&self.state.partial_densities);
        output.state.momentum = update.momentum;
        output.state.energy = update.energy;
        output.volume_m3 = self.volume_m3;
        Ok(update)
    }
}

#[cfg(test)]
mod tests;
