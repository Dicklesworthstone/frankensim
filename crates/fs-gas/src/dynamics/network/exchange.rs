//! Rigid compact-storage candidates from complete extensive port exchanges.
//!
//! The spatial/time owner supplies the opposite mass and TOTAL-energy transfer
//! of its port flux. This leaf retains the original plenum and gas source; it
//! introduces neither a clock nor a second pressure or enthalpy law.

use super::Plenum;
use crate::dynamics::{Conserved, Error, FrozenGas, poll};

impl Plenum {
    /// Prepare a rigid-plenum state after an extensive port exchange.
    ///
    /// Positive `mass_into_kg` [kg] and `energy_into_j` [J] enter this volume.
    /// Energy is the complete transported TOTAL energy in `gas`'s retained
    /// reference, including incoming kinetic energy and flow pressure work.
    /// It becomes internal storage because the plenum has no bulk momentum.
    /// Do not add another enthalpy or pressure-work term. The composing owner
    /// records the opposite momentum transfer at the fixed support separately.
    ///
    /// The caller must use this same frozen source and energy gauge for the
    /// original state and exchange. Both endpoints pass its original caloric
    /// admission. Exhaustion, an unsupported temperature or an unrepresentable
    /// inventory refuses; no density, pressure or energy is clipped.
    ///
    /// Volume remains fixed, zero transfer retains the original state exactly,
    /// and the method allocates no storage. This is an algebraic stage candidate:
    /// the composing owner chooses its time quadrature/CFL and publishes all
    /// coupled owners together. A refusal or cancellation never changes `self`.
    pub fn exchange_candidate(
        &self,
        gas: &FrozenGas,
        mass_into_kg: f64,
        energy_into_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if !mass_into_kg.is_finite() || !energy_into_j.is_finite() {
            return Err(Error::Invalid("finite extensive plenum exchange required"));
        }
        let original = gas.primitive(self.state, &mut cancelled)?;
        let state = Conserved {
            density: if mass_into_kg == 0.0 {
                self.state.density
            } else {
                (self.volume_m3 * self.state.density + mass_into_kg) / self.volume_m3
            },
            momentum: 0.0,
            energy: if energy_into_j == 0.0 {
                self.state.energy
            } else {
                (self.volume_m3 * self.state.energy + energy_into_j) / self.volume_m3
            },
        };
        let candidate = Self::new(self.volume_m3, state)?;
        gas.decode(state, Some(original.temperature), &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(candidate)
    }
}

#[cfg(test)]
mod tests;
