//! Bidirectional ballistic gas exchange on one unchanged network snapshot.
//!
//! Both molecular beams spend their ORIGINAL donor inventories. The original
//! species owner transports bulk partial enthalpy; the separate signed energy
//! correction changes each beam to h-RT/(2M) before ONE endpoint reconstruction.
//! A zero net mass flux can still exchange heat. No pressure-order gate belongs
//! on this temperature-driven path. Closed-network endpoint entropy remains
//! authoritative for the finite proposal.
//!
//! Each interface supplies its explicit thin/free-molecular qualification.
//! Do not add another gas law for the same aperture or use this interface as a
//! long accommodating pore without a separate wall-energy/transport owner.

use std::collections::BTreeMap;

use super::super::super::super::snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS};
use super::super::super::{
    Sum,
    effusion::{EffusiveFlux, EffusiveInterface},
};
use super::super::{SpeciesNetworkHeatTransfer, SpeciesNetworkTransfer};
use crate::thermochemical_gas::energy::{GasEnergyError, parcel::ClosedGasParcel};

#[derive(Clone, Copy, Debug)]
pub struct EffusionNetworkLink<'a> {
    pub id: u64,
    /// Oriented left/right finite reservoirs.
    pub parcels: [u64; 2],
    pub interface: &'a EffusiveInterface,
    /// Total exposed area of separately qualified thin apertures [m2].
    /// Zero closes the connection without querying inactive gas sources.
    pub area_m2: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffusionNetworkLinkStep {
    pub id: u64,
    pub parcels: [u64; 2],
    pub flux: EffusiveFlux,
}

#[derive(Clone, Debug, PartialEq)]
pub enum EffusionNetworkError {
    Invalid(&'static str),
    Face { id: u64, error: GasEnergyError },
    Cancelled,
}
impl core::fmt::Display for EffusionNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-reservoir effusion network: {self:?}")
    }
}
impl std::error::Error for EffusionNetworkError {}
use EffusionNetworkError as Error;

fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn face_error(id: u64, error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled {
        Error::Cancelled
    } else {
        Error::Face { id, error }
    }
}

/// Sample all independent paths without publishing gas, species or energy.
pub(super) fn sample_transfers(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    links: &[EffusionNetworkLink<'_>],
    duration_s: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<
    (
        Vec<SpeciesNetworkTransfer>,
        Vec<SpeciesNetworkHeatTransfer>,
        Vec<EffusionNetworkLinkStep>,
    ),
    Error,
> {
    check(&mut cancelled)?;
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || links.len() > MAX_SNAPSHOT_TRANSFERS
    {
        return Err(Error::Invalid(
            "bounded nonempty effusion graph and positive duration required",
        ));
    }
    let mut ordered = BTreeMap::new();
    for link in links {
        check(&mut cancelled)?;
        if link.parcels[0] == link.parcels[1]
            || link.parcels.iter().any(|id| !parcels.contains_key(id))
            || !link.area_m2.is_finite()
            || link.area_m2 < 0.0
            || ordered.insert(link.id, link).is_some()
        {
            return Err(Error::Invalid(
                "unique effusion IDs, distinct existing endpoints and nonnegative area required",
            ));
        }
    }
    let mut transfers = Vec::new();
    let mut heat = Vec::new();
    let mut reports = Vec::with_capacity(links.len());
    for (id, link) in ordered {
        check(&mut cancelled)?;
        let [left, right] = link.parcels;
        let flux = parcels[&left]
            .effusive_flux_to(
                &parcels[&right],
                link.interface,
                link.area_m2,
                &mut cancelled,
            )
            .map_err(|error| face_error(id, error))?;
        let tl = parcels[&left].state().state().temperature;
        let tr = parcels[&right].state().state().temperature;
        let mut additional = Sum::default();
        for species in &flux.species {
            let mut masses = [0.0; 2];
            for (donor, receiver, side, beam) in [
                (left, right, 0, &species.left_to_right),
                (right, left, 1, &species.right_to_left),
            ] {
                check(&mut cancelled)?;
                let mass = beam.mass_rate_kg_s * duration_s;
                if !mass.is_finite() || mass < 0.0 || (beam.mass_rate_kg_s > 0.0 && mass == 0.0) {
                    return Err(Error::Invalid(
                        "integrated effusive beam mass is unrepresentable",
                    ));
                }
                masses[side] = mass;
                if mass == 0.0 {
                    continue;
                }
                if transfers.len() + heat.len() == MAX_SNAPSHOT_TRANSFERS {
                    return Err(Error::Invalid(
                        "effusive beams exceed the shared transfer bound",
                    ));
                }
                let specific = beam.additional_specific_energy_j_kg.ok_or(Error::Invalid(
                    "active molecular beam lost its energy convention",
                ))?;
                let correction = mass * specific;
                if !correction.is_finite() || correction >= 0.0 {
                    return Err(Error::Invalid(
                        "nonzero effusive molecular energy correction is unrepresentable",
                    ));
                }
                transfers.push(SpeciesNetworkTransfer {
                    donor,
                    receiver,
                    species: species.species.clone(),
                    mass_kg: mass,
                });
            }
            // Apply the correction to the actual rounded masses spent above.
            // Pair opposite beams before taking their temperature difference;
            // subtracting two large products can erase near-equilibrium heat.
            let delta = masses[0] - masses[1];
            if delta != 0.0 {
                let donor = if delta > 0.0 {
                    &species.left_to_right
                } else {
                    &species.right_to_left
                };
                let specific = donor.additional_specific_energy_j_kg.ok_or(Error::Invalid(
                    "net molecular beam lost its energy convention",
                ))?;
                let correction = delta * specific;
                if !correction.is_finite() || correction == 0.0 {
                    return Err(Error::Invalid(
                        "net molecular energy correction is unrepresentable",
                    ));
                }
                additional
                    .add(correction)
                    .map_err(|error| face_error(id, error))?;
            }
            let paired = masses[0].min(masses[1]);
            if paired > 0.0 && tl != tr {
                let specific = species
                    .left_to_right
                    .additional_specific_energy_j_kg
                    .ok_or(Error::Invalid(
                        "paired molecular beam lost its energy convention",
                    ))?;
                let correction = (paired * specific) * ((tl - tr) / tl);
                if !correction.is_finite() || correction == 0.0 {
                    return Err(Error::Invalid(
                        "paired molecular energy correction is unrepresentable",
                    ));
                }
                additional
                    .add(correction)
                    .map_err(|error| face_error(id, error))?;
            }
        }
        if additional.value != 0.0 {
            if transfers.len() + heat.len() == MAX_SNAPSHOT_TRANSFERS {
                return Err(Error::Invalid(
                    "effusive energy corrections exceed the shared transfer bound",
                ));
            }
            heat.push(SpeciesNetworkHeatTransfer {
                parcels: link.parcels,
                heat_left_to_right_j: additional.value,
            });
        }
        reports.push(EffusionNetworkLinkStep {
            id,
            parcels: link.parcels,
            flux,
        });
    }
    check(&mut cancelled)?;
    Ok((transfers, heat, reports))
}

#[cfg(test)]
mod tests;
