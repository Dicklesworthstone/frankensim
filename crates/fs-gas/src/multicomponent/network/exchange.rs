//! Source-admitted rigid mixture storage from extensive port exchanges.
//!
//! Species enter in canonical source order. Complete incoming total energy,
//! including kinetic energy, becomes internal energy at zero bulk momentum.

use super::Plenum;
use crate::multicomponent::{
    Conserved, Error, MAX_COMPONENTS, MixtureGas, allocated, poll, positive,
};

impl Plenum {
    /// Prepare a rigid-plenum candidate from species masses [kg] and energy [J].
    ///
    /// Positive exchanges enter the tank. `species_masses_into_kg` follows
    /// `gas`'s exact canonical species order, and `energy_into_j` uses its
    /// per-species reference gauge. Supply complete transported TOTAL energy;
    /// incoming kinetic energy thermalizes and flow pressure work is already
    /// included. The composing owner records the separate support impulse.
    ///
    /// Both endpoints are admitted by the original source. No clipping,
    /// reaction, moving-volume work, clock or time quadrature is inferred.
    /// Failure leaves `self` unchanged. This convenience form allocates one
    /// bounded species vector; repeated stages can use
    /// [`Self::exchange_candidate_into`] with an existing output plenum.
    pub fn exchange_candidate(
        &self,
        gas: &MixtureGas,
        species_masses_into_kg: &[f64],
        energy_into_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let values =
            self.exchange_values(gas, species_masses_into_kg, energy_into_j, &mut cancelled)?;
        let n = gas.species_count();
        let mut partial_densities = allocated(n, 0.0)?;
        partial_densities.copy_from_slice(&values[..n]);
        let candidate = Self::new(
            self.volume_m3,
            Conserved {
                partial_densities,
                momentum: 0.0,
                energy: values[n + 1],
            },
        )?;
        poll(&mut cancelled)?;
        Ok(candidate)
    }

    /// Write the same candidate into an already allocated plenum.
    ///
    /// The output needs exactly one density slot per source species. Its whole
    /// state and volume are replaced only after all admission and cancellation
    /// checks succeed. Failure preserves every output value and its capacity.
    /// The bounded stack candidate makes successful warm calls allocation-free.
    pub fn exchange_candidate_into(
        &self,
        gas: &MixtureGas,
        species_masses_into_kg: &[f64],
        energy_into_j: f64,
        output: &mut Self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(&mut cancelled)?;
        let n = gas.species_count();
        if output.state.partial_densities.len() != n {
            return Err(Error::Invalid(
                "complete mixture plenum output species axis required",
            ));
        }
        let values =
            self.exchange_values(gas, species_masses_into_kg, energy_into_j, &mut cancelled)?;
        poll(&mut cancelled)?;
        output.state.partial_densities.copy_from_slice(&values[..n]);
        output.state.momentum = 0.0;
        output.state.energy = values[n + 1];
        output.volume_m3 = self.volume_m3;
        Ok(())
    }

    fn exchange_values(
        &self,
        gas: &MixtureGas,
        species_masses_into_kg: &[f64],
        energy_into_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; MAX_COMPONENTS], Error> {
        poll(cancelled)?;
        let n = gas.species_count();
        if species_masses_into_kg.len() != n
            || species_masses_into_kg.iter().any(|x| !x.is_finite())
            || !energy_into_j.is_finite()
        {
            return Err(Error::Invalid(
                "finite canonical-species plenum exchange required",
            ));
        }
        let original = gas.primitive(self.state(), &mut *cancelled)?;
        let mut values = [0.0; MAX_COMPONENTS];
        for i in 0..n {
            poll(cancelled)?;
            let before = self.state.partial_densities[i];
            values[i] = if species_masses_into_kg[i] == 0.0 {
                before
            } else {
                (self.volume_m3 * before + species_masses_into_kg[i]) / self.volume_m3
            };
            if !values[i].is_finite()
                || values[i] < 0.0
                || !(self.volume_m3 * values[i]).is_finite()
            {
                return Err(Error::Invalid(
                    "mixture plenum species inventory exhausted or unrepresentable",
                ));
            }
        }
        values[n + 1] = if energy_into_j == 0.0 {
            self.state.energy
        } else {
            (self.volume_m3 * self.state.energy + energy_into_j) / self.volume_m3
        };
        if !(self.volume_m3 * values[n + 1]).is_finite() {
            return Err(Error::Invalid(
                "unrepresentable mixture plenum energy inventory",
            ));
        }
        let admitted = gas.decode(&values[..n + 2], Some(original.temperature), cancelled)?;
        if !positive(self.volume_m3 * admitted.density) {
            return Err(Error::Invalid("unrepresentable total mixture plenum mass"));
        }
        poll(cancelled)?;
        Ok(values)
    }
}

#[cfg(test)]
mod tests;
