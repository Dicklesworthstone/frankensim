//! Finite gases exchanging through molecular channels and held thermal walls.
//!
//! Every transmitted and reflected stream uses one original gas snapshot and
//! shares the existing gross-withdrawal allowance. Accommodation heat and the
//! original molecular crossing energy enter ONE caloric solve per gas. The base
//! step uses prescribed thermostats with explicit opposite heat and entropy.
//! Gas plus thermostat energy/entropy is audited before publication. The thermal
//! submodule instead owns finite material-wall storage and actual entropy.
//! This is a first-order, quasistatic channel stage; callers refine time.

/// Finite material walls sharing the original fixed-volume caloric owner.
pub mod thermal;

use super::{EffusionNetworkLink, check, face_error};
use crate::thermochemical_gas::energy::parcel::transfer::{
    snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
    species::{
        diffusion::network::DiffusionNetworkError,
        effusion::{
            EffusiveFlux,
            channel::{ChannelWallExchange, DiffuseMolecularChannel},
        },
        entropy,
        network::{
            SpeciesNetworkExternalHeat, SpeciesNetworkStep,
            boundary::{BoundarySpeciesTransfer, SpeciesNetworkBoundary},
            exchange_candidate_with_external_heat,
            transport::{
                self, GasTransportError, GasTransportPaths, GasTransportRates, GasTransportSolve,
            },
        },
    },
};
use crate::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_math::dd::Dd;
use std::collections::BTreeMap;

/// Separate circular channels, with a common held temperature along each bore.
#[derive(Clone, Copy, Debug)]
pub struct DiffuseChannelNetworkLink<'a> {
    pub id: u64,
    pub parcels: [u64; 2],
    pub channel: &'a DiffuseMolecularChannel,
    /// Sum of physical mouth areas, independently of individual bore dimensions.
    pub area_m2: f64,
    pub wall_temperature_k: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DiffuseChannelLinkStep {
    pub id: u64,
    pub parcels: [u64; 2],
    pub channel: DiffuseMolecularChannel,
    pub wall_temperature_k: f64,
    pub transmitted: EffusiveFlux,
    pub wall: ChannelWallExchange,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DiffuseChannelNetworkStep {
    pub duration_s: f64,
    pub transport: GasTransportRates,
    pub channels: Vec<DiffuseChannelLinkStep>,
    pub exchange: SpeciesNetworkStep,
    pub heat_into_walls_j: f64,
    /// Independent original-to-final gas storage plus opposite thermostat heat.
    pub energy_residual_j: f64,
    pub gas_entropy_change_j_k: f64,
    pub wall_entropy_change_j_k: f64,
    pub total_entropy_change_j_k: f64,
}
impl DiffuseChannelNetworkStep {
    /// Recheck original gas/collision sources and held-wall caloric domains.
    /// Newly arriving species must admit the actual wall before publication.
    pub fn admit_endpoint_driving_forces(
        &self,
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), GasTransportError> {
        self.transport
            .admit_endpoint_driving_forces(parcels, &mut cancelled)?;
        for row in &self.channels {
            check(&mut cancelled)?;
            if row.transmitted.species.is_empty() {
                continue;
            }
            for id in row.parcels {
                let parcel = parcels.get(&id).ok_or(GasTransportError::Invalid(
                    "channel endpoint lost an original gas identity",
                ))?;
                row.channel
                    .admit_wall_for_parcel(parcel, row.wall_temperature_k, &mut cancelled)
                    .map_err(|error| face_error(row.id, error))?;
            }
        }
        check(&mut cancelled)?;
        Ok(())
    }
}

/// Complete candidate; caller-owned gas inventories never change on refusal.
#[derive(Clone, Debug, PartialEq)]
pub struct DiffuseChannelNetworkCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: DiffuseChannelNetworkStep,
}
impl DiffuseChannelNetworkCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }
    pub const fn step(&self) -> &DiffuseChannelNetworkStep {
        &self.step
    }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, DiffuseChannelNetworkStep) {
        (self.parcels, self.step)
    }
}

/// Advance native transport and fully accommodating channels simultaneously.
///
/// A held wall supplies/absorbs heat; this is not a thermally isolated gas pair.
/// Reflected molecules are represented as equal incoming/outgoing species rows
/// with zero net mass/enthalpy, while their accommodation heat is counted once.
/// This spends the original incident inventory together with all other paths.
/// No continuum nozzle pressure work or second mixing/thermalization heat enters.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    paths: GasTransportPaths<'_>,
    links: &[DiffuseChannelNetworkLink<'_>],
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: GasTransportSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DiffuseChannelNetworkCandidate, GasTransportError> {
    check(&mut cancelled)?;
    let count = [
        paths.nozzles.len(),
        paths.diffusion.len(),
        paths.porous.len(),
        paths.membranes.len(),
        paths.dusty_gas.len(),
        paths.layered_dusty_gas.len(),
        paths.effusion.len(),
        links.len(),
    ]
    .into_iter()
    .try_fold(0usize, |sum, n| sum.checked_add(n));
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || count.is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS)
        || !controls.entropy_tolerance_j_k.is_finite()
        || controls.entropy_tolerance_j_k <= 0.0
    {
        return Err(GasTransportError::Invalid(
            "channel transport needs bounded gases/paths, positive duration and entropy allowance",
        ));
    }
    let mut ordered = BTreeMap::new();
    for link in links {
        check(&mut cancelled)?;
        if link.parcels[0] == link.parcels[1]
            || link.parcels.iter().any(|id| !parcels.contains_key(id))
            || ordered.insert(link.id, link).is_some()
        {
            return Err(GasTransportError::Invalid(
                "channels need unique IDs and two distinct existing gases",
            ));
        }
    }
    let mut molecular = Vec::with_capacity(links.len());
    for link in ordered.values() {
        check(&mut cancelled)?;
        let area_m2 = link
            .channel
            .transmitted_area_m2(link.area_m2)
            .map_err(|error| face_error(link.id, error))?;
        if area_m2 > 0.0 {
            for id in link.parcels {
                link.channel
                    .admit_parcel(&parcels[&id], &mut cancelled)
                    .map_err(|error| face_error(link.id, error))?;
            }
        }
        molecular.push(EffusionNetworkLink {
            id: link.id,
            parcels: link.parcels,
            interface: link.channel.mouth(),
            area_m2,
        });
    }
    let mut sampled = transport::sample_transfers(
        parcels,
        paths,
        duration_s,
        endpoints,
        controls,
        &mut cancelled,
    )?;
    let (transfers, heat, reports) =
        super::sample_transfers(parcels, &molecular, duration_s, &mut cancelled)?;
    if sampled.transfers.len() + sampled.heat.len() + transfers.len() + heat.len()
        > MAX_SNAPSHOT_TRANSFERS
    {
        return Err(GasTransportError::Invalid(
            "channel streams exceed the shared work bound",
        ));
    }
    sampled.transfers.extend(transfers);
    sampled.heat.extend(heat);
    let mut channels = Vec::with_capacity(links.len());
    let mut returned = Vec::with_capacity(2 * links.len());
    let mut external_heat = Vec::with_capacity(2 * links.len());
    let mut wall_heat = Dd::ZERO;
    let mut wall_entropy = Dd::ZERO;
    let mut work = sampled.transfers.len() + sampled.heat.len();
    for (index, (link, report)) in ordered.values().zip(reports).enumerate() {
        check(&mut cancelled)?;
        let wall = link
            .channel
            .integrate_wall_exchange(
                link.parcels.map(|id| &parcels[&id]),
                &report.flux,
                link.wall_temperature_k,
                duration_s,
                &mut cancelled,
            )
            .map_err(|error| face_error(link.id, error))?;
        for side in 0..2 {
            let id = (2 * index + side) as u64;
            let mut rows = Vec::with_capacity(wall.reflected_species.len());
            for reflected in &wall.reflected_species {
                check(&mut cancelled)?;
                let mass = reflected.mass_kg[side];
                if mass > 0.0 {
                    rows.push(BoundarySpeciesTransfer {
                        species: reflected.species.clone(),
                        incoming_mass_kg: mass,
                        outgoing_mass_kg: mass,
                        net_mass_into_gas_kg: 0.0,
                        relative_energy_into_gas_j: 0.0,
                    });
                }
            }
            if !rows.is_empty() {
                work += rows.len() + 1;
                returned.push((
                    id,
                    link.parcels[side],
                    link.channel.mouth().gas_catalog(),
                    rows,
                ));
            }
            if wall.heat_into_gas_j[side] != 0.0 {
                work += 1;
                external_heat.push(SpeciesNetworkExternalHeat {
                    id,
                    parcel: link.parcels[side],
                    heat_into_gas_j: wall.heat_into_gas_j[side],
                });
            }
        }
        if work > MAX_SNAPSHOT_TRANSFERS {
            return Err(GasTransportError::Invalid(
                "channel wall/return streams exceed the shared work bound",
            ));
        }
        wall_heat = wall_heat + Dd::from_f64(wall.heat_into_wall_j);
        if wall.heat_into_wall_j != 0.0 {
            wall_entropy =
                wall_entropy + Dd::from_f64(wall.heat_into_wall_j / link.wall_temperature_k);
        }
        channels.push(DiffuseChannelLinkStep {
            id: link.id,
            parcels: link.parcels,
            channel: link.channel.clone(),
            wall_temperature_k: link.wall_temperature_k,
            transmitted: report.flux,
            wall,
        });
    }
    let boundaries: Vec<_> = returned
        .iter()
        .map(|(id, parcel, catalog, rows)| SpeciesNetworkBoundary {
            id: *id,
            parcel: *parcel,
            catalog,
            species: rows,
        })
        .collect();
    let candidate = exchange_candidate_with_external_heat(
        parcels,
        &sampled.transfers,
        &sampled.heat,
        &boundaries,
        &external_heat,
        endpoints,
        controls.exchange,
        &mut cancelled,
    )?;
    // Actual storage and gas entropy, with original formation datums removed
    // before arithmetic. Thermostat heat is outside the closed-gas audit itself.
    let change = entropy::closed_inventory_change(
        parcels
            .iter()
            .map(|(&id, old)| (id, old, &candidate.parcels()[&id])),
        controls.exchange.mass_tolerance_kg,
        &mut cancelled,
    )
    .map_err(|error| match error {
        entropy::Error::Invalid(reason) => GasTransportError::Invalid(reason),
        entropy::Error::Parcel { id, error } => {
            GasTransportError::Diffusion(DiffusionNetworkError::Parcel { id, error })
        }
        entropy::Error::Cancelled => GasTransportError::Cancelled,
    })?;
    let energy_residual_j = (Dd::from_f64(change.energy_change_j) + wall_heat).to_f64();
    let total_entropy_change_j_k =
        (Dd::from_f64(change.entropy_change_j_k) + wall_entropy).to_f64();
    if !energy_residual_j.is_finite()
        || !total_entropy_change_j_k.is_finite()
        || change.total_mass_residual_kg.abs() > controls.exchange.mass_tolerance_kg
        || change.maximum_species_mass_residual_kg > controls.exchange.mass_tolerance_kg
        || energy_residual_j.abs() > controls.exchange.energy_tolerance_j
    {
        return Err(GasTransportError::Invalid(
            "channel gas/wall conservation exceeded its fixed budget",
        ));
    }
    if total_entropy_change_j_k < -controls.entropy_tolerance_j_k {
        return Err(GasTransportError::Diffusion(
            DiffusionNetworkError::EntropyDecrease {
                change_j_k: total_entropy_change_j_k,
                tolerance_j_k: controls.entropy_tolerance_j_k,
            },
        ));
    }
    let (parcels, exchange) = candidate.into_parts();
    let step = DiffuseChannelNetworkStep {
        duration_s,
        transport: sampled.rates,
        channels,
        exchange,
        heat_into_walls_j: wall_heat.to_f64(),
        energy_residual_j,
        gas_entropy_change_j_k: change.entropy_change_j_k,
        wall_entropy_change_j_k: wall_entropy.to_f64(),
        total_entropy_change_j_k,
    };
    step.admit_endpoint_driving_forces(&parcels, &mut cancelled)?;
    check(&mut cancelled)?;
    Ok(DiffuseChannelNetworkCandidate { parcels, step })
}

#[cfg(test)]
mod tests;
