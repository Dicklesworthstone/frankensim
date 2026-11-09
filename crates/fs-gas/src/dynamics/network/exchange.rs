//! Rigid compact-storage candidates from complete extensive port exchanges.
//!
//! The spatial/time owner supplies the opposite mass and TOTAL-energy transfer
//! of its port flux. This leaf retains the original plenum and gas source; it
//! introduces neither a clock nor a second pressure or enthalpy law.

use super::{Plenum, PlenumBoundaryWork};
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

    /// Prepare extensive port exchange and prescribed chamber-volume work.
    ///
    /// `final_volume_m3` is the represented physical endpoint [m³]. Port mass
    /// [kg] and TOTAL energy [J] have the same sign and source gauge as
    /// [`Self::exchange_candidate`]. Moving-wall work is additional to that
    /// transported energy; the port already contains its own flow work.
    ///
    /// With no `forward_stage`, the source pressure at `self` supplies a
    /// forward-Euler candidate. To finish a native SSP RK2 interval, supply
    /// that complete forward candidate and the two half-weighted port
    /// transfers. The work is then `-0.5 * (p_initial + p_forward) * delta_V`.
    /// Both pressures come from this original gas, and the forward stage must
    /// have the same endpoint volume. The composing owner supplies one affine
    /// volume interval, enforces its acoustic/volume CFL and publishes the
    /// field, chamber and mechanical receipt together.
    ///
    /// Extensive mass and internal energy are divided by the final volume;
    /// densities at different volumes are never averaged. Work uses the
    /// represented volume difference, so an unrepresented stroke cannot spend
    /// energy. Equal volumes take the original rigid path exactly. No clock,
    /// heat source, piston inertia or extra caloric law is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn exchange_with_volume_candidate(
        &self,
        gas: &FrozenGas,
        final_volume_m3: f64,
        mass_into_kg: f64,
        energy_into_j: f64,
        forward_stage: Option<&Self>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, PlenumBoundaryWork), Error> {
        poll(&mut cancelled)?;
        if !final_volume_m3.is_finite()
            || final_volume_m3 <= 0.0
            || !mass_into_kg.is_finite()
            || !energy_into_j.is_finite()
            || forward_stage.is_some_and(|stage| stage.volume_m3 != final_volume_m3)
        {
            return Err(Error::Invalid(
                "moving plenum needs finite transfers and a common positive endpoint volume",
            ));
        }
        let mut work = PlenumBoundaryWork {
            initial_volume_m3: self.volume_m3,
            final_volume_m3,
            work_on_gas_j: 0.0,
        };
        if final_volume_m3 == self.volume_m3 {
            return Ok((
                self.exchange_candidate(gas, mass_into_kg, energy_into_j, cancelled)?,
                work,
            ));
        }
        let initial = gas.primitive(self.state, &mut cancelled)?;
        let pressure = if let Some(stage) = forward_stage {
            let forward = gas.primitive(stage.state, &mut cancelled)?;
            0.5 * initial.pressure + 0.5 * forward.pressure
        } else {
            initial.pressure
        };
        work.work_on_gas_j = -pressure * (final_volume_m3 - self.volume_m3);
        if !work.work_on_gas_j.is_finite() || work.work_on_gas_j == 0.0 {
            return Err(Error::Invalid(
                "unrepresentable moving-plenum pressure work",
            ));
        }
        let state = Conserved {
            density: (self.volume_m3 * self.state.density + mass_into_kg) / final_volume_m3,
            momentum: 0.0,
            energy: (self.volume_m3 * self.state.energy + energy_into_j + work.work_on_gas_j)
                / final_volume_m3,
        };
        let candidate = Self::new(final_volume_m3, state)?;
        gas.decode(state, Some(initial.temperature), &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok((candidate, work))
    }
}

#[cfg(test)]
mod tests;
