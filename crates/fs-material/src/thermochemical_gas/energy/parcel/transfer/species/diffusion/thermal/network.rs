//! Thermal diffusion joins all original gas paths before ONE inventory update.
//!
//! Each link has reusable constitutive scratch. Every rate samples the original
//! state; donor species enthalpy and the source-conjugate heat correction enter
//! the existing conservative owner together. Incoming material cannot finance
//! an outgoing allowance. Final source support and total network entropy are
//! admitted before any endpoint is returned. Caller-controlled time refinement
//! remains necessary; this explicit update is not an implicit stability claim.

use super::{
    ClosedGasParcel, GasDiffusionGeometry, GasEnergyError, GasThermodiffusionFlux,
    GasThermodiffusionModel, GasThermodiffusionSolve, GasThermodiffusionWorkspace, positive,
};
use crate::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::transfer::{
        snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
        species::network::{
            SpeciesNetworkHeatTransfer, SpeciesNetworkTransfer,
            transport::{
                self, GasTransportCandidate, GasTransportError, GasTransportPaths,
                GasTransportSolve,
            },
        },
    },
};
use std::collections::{BTreeMap, BTreeSet};

/// An explicitly exposed molecular face, positive from parcels[0] to parcels[1].
/// It already includes ordinary molecular diffusion and heat conduction; do not
/// add either mechanism again for this same physical exposure in `other_paths`.
#[derive(Clone, Debug)]
pub struct GasThermodiffusionLink<'a> {
    pub id: u64,
    pub parcels: [u64; 2],
    pub geometry: GasDiffusionGeometry,
    pub model: &'a GasThermodiffusionModel<'a>,
    pub solve: GasThermodiffusionSolve,
}
#[derive(Clone, Debug, PartialEq)]
pub struct GasThermodiffusionLinkStep {
    pub id: u64,
    pub parcels: [u64; 2],
    pub flux: GasThermodiffusionFlux,
}
#[derive(Clone, Debug, PartialEq)]
pub enum GasThermodiffusionNetworkError {
    Invalid(&'static str),
    Face { id: u64, error: GasEnergyError },
    Transport(GasTransportError),
    Cancelled,
}
impl core::fmt::Display for GasThermodiffusionNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "gas thermodiffusion network: {self:?}")
    }
}
impl std::error::Error for GasThermodiffusionNetworkError {}

/// Source-bearing thermal faces plus the existing complete conservative result.
#[derive(Clone, Debug, PartialEq)]
pub struct GasThermodiffusionNetworkCandidate {
    transport: GasTransportCandidate,
    links: Vec<GasThermodiffusionLinkStep>,
}
impl GasThermodiffusionNetworkCandidate {
    pub const fn transport(&self) -> &GasTransportCandidate {
        &self.transport
    }
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        self.transport.parcels()
    }
    pub fn links(&self) -> &[GasThermodiffusionLinkStep] {
        &self.links
    }
    pub fn into_parts(self) -> (GasTransportCandidate, Vec<GasThermodiffusionLinkStep>) {
        (self.transport, self.links)
    }
}
use GasThermodiffusionNetworkError as Error;
fn face_error(id: u64, error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled {
        Error::Cancelled
    } else {
        Error::Face { id, error }
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

/// Advance nonisothermal molecular faces together with any other ORIGINAL paths.
/// Workspaces correspond to input links and may be reused across all timesteps.
/// Link IDs order deterministic sampling, never sequential physical updates.
/// A failure changes only scratch; all parcel inputs remain unchanged.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    links: &[GasThermodiffusionLink<'_>],
    workspaces: &mut [GasThermodiffusionWorkspace],
    other_paths: GasTransportPaths<'_>,
    duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: GasTransportSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GasThermodiffusionNetworkCandidate, Error> {
    poll(&mut cancelled)?;
    if !positive(duration_s)
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || links.len() != workspaces.len()
        || links.len() > MAX_SNAPSHOT_TRANSFERS
        || links
            .iter()
            .try_fold(0usize, |n, link| {
                n.checked_add(link.model.species().len() + 1)
            })
            .is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS)
    {
        return Err(Error::Invalid(
            "thermodiffusion snapshot dimensions or work bound",
        ));
    }
    let mut ids = BTreeSet::new();
    for link in links {
        if !ids.insert(link.id)
            || link.parcels[0] == link.parcels[1]
            || link.parcels.iter().any(|id| !parcels.contains_key(id))
        {
            return Err(Error::Invalid(
                "thermodiffusion link identities or endpoints",
            ));
        }
    }
    let mut order: Vec<_> = (0..links.len()).collect();
    order.sort_by_key(|&i| links[i].id);
    let mut transfers = Vec::new();
    let mut heat = Vec::with_capacity(links.len());
    let mut sampled = Vec::with_capacity(links.len());
    for &index in &order {
        poll(&mut cancelled)?;
        let link = &links[index];
        let flux = link
            .model
            .face_flux(
                &parcels[&link.parcels[0]],
                &parcels[&link.parcels[1]],
                link.geometry,
                link.solve,
                &mut workspaces[index],
                &mut cancelled,
            )
            .map_err(|error| face_error(link.id, error))?;
        for row in &flux.species {
            let signed = row.mass_rate_left_to_right_kg_s * duration_s;
            if !signed.is_finite() || (row.mass_rate_left_to_right_kg_s != 0.0 && signed == 0.0) {
                return Err(Error::Invalid(
                    "unrepresentable integrated thermodiffusion mass",
                ));
            }
            if signed != 0.0 {
                let [donor, receiver] = if signed > 0.0 {
                    link.parcels
                } else {
                    [link.parcels[1], link.parcels[0]]
                };
                transfers.push(SpeciesNetworkTransfer {
                    donor,
                    receiver,
                    species: row.species.clone(),
                    mass_kg: signed.abs(),
                });
            }
        }
        let energy = flux.additional_heat_rate_left_to_right_w * duration_s;
        if !energy.is_finite()
            || (flux.additional_heat_rate_left_to_right_w != 0.0 && energy == 0.0)
        {
            return Err(Error::Invalid(
                "unrepresentable integrated thermodiffusion heat",
            ));
        }
        if energy != 0.0 {
            heat.push(SpeciesNetworkHeatTransfer {
                parcels: link.parcels,
                heat_left_to_right_j: energy,
            });
        }
        sampled.push(GasThermodiffusionLinkStep {
            id: link.id,
            parcels: link.parcels,
            flux,
        });
    }
    let transport = transport::step_candidate_with_paths_and_fluxes(
        parcels,
        other_paths,
        duration_s,
        &transfers,
        &heat,
        endpoints,
        controls,
        &mut cancelled,
    )
    .map_err(|error| {
        if error == GasTransportError::Cancelled {
            Error::Cancelled
        } else {
            Error::Transport(error)
        }
    })?;
    for index in order {
        let link = &links[index];
        link.model
            .admit_endpoints(
                &transport.parcels()[&link.parcels[0]],
                &transport.parcels()[&link.parcels[1]],
                link.solve,
                &mut cancelled,
            )
            .map_err(|error| face_error(link.id, error))?;
    }
    poll(&mut cancelled)?;
    Ok(GasThermodiffusionNetworkCandidate {
        transport,
        links: sampled,
    })
}
