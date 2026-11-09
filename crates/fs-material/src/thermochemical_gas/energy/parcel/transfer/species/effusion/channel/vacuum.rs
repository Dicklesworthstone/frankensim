//! Molecular escape through a diffuse channel into an ideal absorbing vacuum.
//!
//! The vacuum supplies no incident beam or thermodynamic state. Of each source
//! incident beam, `tau0` escapes directly, `tau-tau0` escapes after wall
//! scattering and `1-tau` returns to the SAME original reservoir. Full thermal
//! accommodation gives each wall-hit molecule energy change
//! `e_i(Twall)-e_i(Tgas)`, where `e_i(T)=h_i(T)-R*T/(2*M_i)` [J/kg]. The optional
//! partial-accommodation law instead weights this change by the resolved
//! probability of thermalization over ALL collisions before that exit.
//!
//! Returned-molecule heat enters the finite gas. The heat added to escaping
//! molecules belongs to the outward energy ledger, not the gas caloric solve.
//! Opposite wall heat closes both exchanges without inventing wall storage.
//! This inherits the parent channel's collisionless, qualified accommodation,
//! stationary-wall assumptions and estimated geometric transmission. No pump
//! work, vacuum temperature, receiver inventory or closed-system entropy claim
//! is supplied.

use super::super::vacuum::VacuumEffusiveFlux;
use super::{
    ClosedGasParcel, DiffuseMolecularChannel, GasEnergyError, SpeciesId, invalid, poll, positive,
    positive_product_ratio,
};
use fs_math::dd::Dd;

/// Gross molecular mass that leaves and returns to the same original gas.
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumChannelReflectedSpecies {
    pub species: SpeciesId,
    pub mass_kg: f64,
}

/// Integrated accommodation in addition to original transmitted molecular energy.
///
/// Positive heat enters the finite gas, escaped molecules or wall thermostat.
/// The gas loses only the original transmitted molecular energy and receives
/// `heat_into_gas_j`; `heat_into_escaped_j` is supplied by the wall after escape
/// from the source gas. The caller owns the thermostat's opposite heat.
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumChannelWallExchange {
    pub reflected_species: Vec<VacuumChannelReflectedSpecies>,
    pub heat_into_gas_j: f64,
    pub heat_into_escaped_j: f64,
    pub heat_into_wall_j: f64,
}

/// Gross return rate; it spends the source's original withdrawal allowance.
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumChannelReflectedSpeciesRate {
    pub species: SpeciesId,
    pub mass_rate_kg_s: f64,
}

/// Wall power beyond the original transmitted molecular energy rate [W].
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumChannelWallRates {
    pub reflected_species: Vec<VacuumChannelReflectedSpeciesRate>,
    pub heat_into_gas_w: f64,
    pub heat_into_escaped_w: f64,
    pub heat_into_wall_w: f64,
}

/// One-way rates from a qualified finite gas through a channel into vacuum.
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumDiffuseChannelFlux {
    /// Original aperture law at the channel's transmitted area. Its species
    /// masses are the escaped masses; its energies retain the SOURCE gas datum
    /// and temperature before the separate wall accommodation below.
    pub transmitted: VacuumEffusiveFlux,
    pub wall: VacuumChannelWallRates,
    /// Complete escaped power, including heat exchanged with the wall [W].
    /// Arbitrary formation references can make this negative or obscure a small
    /// wall correction; `wall.heat_into_escaped_w` retains that correction.
    pub escaped_energy_rate_w: f64,
}

impl DiffuseMolecularChannel {
    /// Sample a true-vacuum outlet without constructing a receiving gas.
    ///
    /// Gas power is `-transmitted.energy_rate_w + wall.heat_into_gas_w`.
    /// Escaped power is `escaped_energy_rate_w`. Together with
    /// `wall.heat_into_wall_w` these close the channel energy balance.
    /// Source species and collision qualifications remain active for any
    /// positive area. Zero area and zero length do not query an inactive wall.
    ///
    /// Finite transactions must include BOTH escaped and reflected species in
    /// gross source withdrawals, and integrate the shared accommodation kernel
    /// using the same rounded transmitted masses as the boundary owner.
    pub fn sample_to_vacuum(
        &self,
        gas: &ClosedGasParcel,
        area_m2: f64,
        wall_temperature_k: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VacuumDiffuseChannelFlux, GasEnergyError> {
        poll(&mut cancelled)?;
        let area = self.transmitted_area_m2(area_m2)?;
        if area > 0.0 {
            self.admit_parcel(gas, &mut cancelled)?;
        }
        let transmitted = gas.effusive_flux_to_vacuum(self.mouth(), area, &mut cancelled)?;
        let exchange = self.integrate_vacuum_wall_exchange(
            gas,
            &transmitted,
            wall_temperature_k,
            1.0,
            &mut cancelled,
        )?;
        let escaped_energy_rate_w = (Dd::from_f64(transmitted.energy_rate_w)
            + Dd::from_f64(exchange.heat_into_escaped_j))
        .to_f64();
        if !escaped_energy_rate_w.is_finite() {
            return Err(invalid("channel escaped energy rate is unrepresentable"));
        }
        let wall = VacuumChannelWallRates {
            reflected_species: exchange
                .reflected_species
                .into_iter()
                .map(|row| VacuumChannelReflectedSpeciesRate {
                    species: row.species,
                    mass_rate_kg_s: row.mass_kg,
                })
                .collect(),
            heat_into_gas_w: exchange.heat_into_gas_j,
            heat_into_escaped_w: exchange.heat_into_escaped_j,
            heat_into_wall_w: exchange.heat_into_wall_j,
        };
        poll(&mut cancelled)?;
        Ok(VacuumDiffuseChannelFlux {
            transmitted,
            wall,
            escaped_energy_rate_w,
        })
    }

    /// Reuse the ORIGINAL transmitted vacuum beam and native boundary mass
    /// rounding. No second source sample or receiving state is introduced.
    pub(crate) fn integrate_vacuum_wall_exchange(
        &self,
        gas: &ClosedGasParcel,
        transmitted: &VacuumEffusiveFlux,
        wall_temperature_k: f64,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VacuumChannelWallExchange, GasEnergyError> {
        poll(&mut cancelled)?;
        if !positive(duration_s) || transmitted.interface != *self.mouth() {
            return Err(invalid(
                "vacuum channel exchange lost its original interface or positive interval",
            ));
        }
        let mut result = VacuumChannelWallExchange {
            reflected_species: Vec::with_capacity(transmitted.species.len()),
            heat_into_gas_j: 0.0,
            heat_into_escaped_j: 0.0,
            heat_into_wall_j: 0.0,
        };
        if transmitted.species.is_empty() || self.geometry().length_m() == 0.0 {
            poll(&mut cancelled)?;
            return Ok(result);
        }
        if self.wall_heat_active() {
            self.admit_wall_temperature(wall_temperature_k)?;
        }
        let temperature = gas.state().state().temperature;
        let (mut into_gas, mut into_escaped) = (Dd::ZERO, Dd::ZERO);
        for row in &transmitted.species {
            poll(&mut cancelled)?;
            if row.outward.mass_rate_kg_s == 0.0 {
                continue;
            }
            // The native vacuum owner uses this exact normalized product for
            // the mass in its boundary row. Energy must use that SAME mass.
            let mass = positive_product_ratio(&[row.outward.mass_rate_kg_s, duration_s], &[])?;
            let wall = self.integrate_beam_wall_exchange(
                self.original_species(&row.species)?,
                mass,
                temperature,
                wall_temperature_k,
            )?;
            if wall.reflected_mass_kg > 0.0 {
                result
                    .reflected_species
                    .push(VacuumChannelReflectedSpecies {
                        species: row.species.clone(),
                        mass_kg: wall.reflected_mass_kg,
                    });
            }
            into_gas = into_gas + Dd::from_f64(wall.heat_into_source_j);
            into_escaped = into_escaped + Dd::from_f64(wall.heat_into_transmitted_j);
        }
        result.heat_into_gas_j = into_gas.to_f64();
        result.heat_into_escaped_j = into_escaped.to_f64();
        result.heat_into_wall_j = (-(into_gas + into_escaped)).to_f64();
        if [
            result.heat_into_gas_j,
            result.heat_into_escaped_j,
            result.heat_into_wall_j,
        ]
        .into_iter()
        .any(|value| !value.is_finite())
        {
            return Err(invalid("vacuum channel wall heat is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
