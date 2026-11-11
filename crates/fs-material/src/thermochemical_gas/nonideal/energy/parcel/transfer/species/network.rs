//! Simultaneous selective streams and heat on original real-gas inventories.
//!
//! Every edge uses its ORIGINAL donor's partial enthalpy. Gross withdrawals
//! are admitted before source mixing, and every affected endpoint is solved
//! once. Incoming material cannot leave again during the same explicit step.
//! Prescribed streams/heat do not supply rates, separation work or an entropy
//! guarantee. A physical transport owner chooses and refines this interval.

pub mod boundary;
mod caloric;
mod snapshot;

use super::{GasEnergyError, GasEnergySolve, IsochoricRealGasParcel, Sum};
use std::collections::BTreeMap;

pub use crate::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkError, SpeciesNetworkExternalHeat, SpeciesNetworkFlux,
    SpeciesNetworkHeatTransfer, SpeciesNetworkParcelStep, SpeciesNetworkSolve,
    SpeciesNetworkTransfer,
};
use SpeciesNetworkError as Error;
use boundary::{RealGasSpeciesNetworkBoundary, SpeciesNetworkBoundaryFlux};

/// Complete original-stream and additional-heat accounts, in canonical order.
/// External heat's opposite finite storage belongs to the enclosing operator.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasSpeciesNetworkStep {
    pub fluxes: Vec<SpeciesNetworkFlux>,
    pub heat: Vec<SpeciesNetworkHeatTransfer>,
    pub external_heat: Vec<SpeciesNetworkExternalHeat>,
    /// External phase streams on the explicitly declared ideal-species
    /// reference chart. Their complete energy is included once in each node.
    pub boundaries: Vec<SpeciesNetworkBoundaryFlux>,
    pub parcels: Vec<SpeciesNetworkParcelStep>,
    pub total_mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    /// Sum of local original-caloric defects, after the species reference
    /// cancels. Includes external heat once, separately from stream enthalpy.
    pub total_energy_residual_j: f64,
}

/// All proposed endpoints, or none. Accepted inputs are never mutated.
#[derive(Clone, Debug, PartialEq)]
pub struct RealGasSpeciesNetworkCandidate {
    parcels: BTreeMap<u64, IsochoricRealGasParcel>,
    step: RealGasSpeciesNetworkStep,
}
impl RealGasSpeciesNetworkCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, IsochoricRealGasParcel> {
        &self.parcels
    }
    pub const fn step(&self) -> &RealGasSpeciesNetworkStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        BTreeMap<u64, IsochoricRealGasParcel>,
        RealGasSpeciesNetworkStep,
    ) {
        (self.parcels, self.step)
    }
}

type Parcels<'a> = BTreeMap<u64, &'a IsochoricRealGasParcel>;
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn material(id: u64, error: GasEnergyError) -> Error {
    if error == GasEnergyError::Cancelled {
        Error::Cancelled
    } else {
        Error::Parcel { id, error }
    }
}

/// Combine all species edges from one unchanged snapshot. Parallel rows are
/// aggregated in deterministic donor/receiver/species/mass order. Source joins
/// contain only the original node and its ORIGINAL incoming donors, so an
/// unrelated disconnected material cannot change another component's domain.
/// Idle nodes retain their exact accepted state. The map and all row kinds
/// share the existing 256-node/4096-submitted-record snapshot bounds.
pub fn exchange_candidate(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
    transfers: &[SpeciesNetworkTransfer],
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<RealGasSpeciesNetworkCandidate, Error> {
    exchange_candidate_with_external_heat(
        parcels,
        transfers,
        &[],
        &[],
        endpoints,
        controls,
        cancelled,
    )
}

/// Additional internal heat shares the same conserved endpoint solve as all
/// species streams. Positive heat leaves parcels[0] and enters parcels[1].
/// `mass_allowance_kg` bounds, per parcel, the rounding of a net mass that the
/// streams' own law sets to zero (a mass-frame diffusion face): a parcel whose
/// computed net mass lies within its allowance keeps its mass exactly. A parcel
/// absent from the map has none, the exact balance every other entry requires.
pub fn exchange_candidate_with_heat(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
    transfers: &[SpeciesNetworkTransfer],
    heat: &[SpeciesNetworkHeatTransfer],
    mass_allowance_kg: &BTreeMap<u64, f64>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<RealGasSpeciesNetworkCandidate, Error> {
    allowed_candidate(
        parcels,
        transfers,
        heat,
        &[],
        &[],
        mass_allowance_kg,
        endpoints,
        controls,
        cancelled,
    )
}

/// Joint species, antisymmetric gas/gas heat, and explicitly identified external
/// heat. External IDs are unique, including zero rows. No heat is hidden in
/// species enthalpy or applied to an already-reconstructed candidate. This
/// supplies no external thermostat, boundary material or entropy audit.
#[allow(clippy::too_many_arguments)]
pub fn exchange_candidate_with_external_heat(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
    transfers: &[SpeciesNetworkTransfer],
    heat: &[SpeciesNetworkHeatTransfer],
    external_heat: &[SpeciesNetworkExternalHeat],
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<RealGasSpeciesNetworkCandidate, Error> {
    exchange_candidate_with_boundaries_and_external_heat(
        parcels,
        transfers,
        heat,
        &[],
        external_heat,
        endpoints,
        controls,
        cancelled,
    )
}

/// External species/phase streams share the original snapshot and its gross
/// withdrawal limits with every internal edge. Boundary energy is explicitly
/// authored on the retained NASA species reference, not inferred latent heat.
#[allow(clippy::too_many_arguments)]
pub fn exchange_candidate_with_boundaries(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
    transfers: &[SpeciesNetworkTransfer],
    heat: &[SpeciesNetworkHeatTransfer],
    boundaries: &[RealGasSpeciesNetworkBoundary<'_>],
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<RealGasSpeciesNetworkCandidate, Error> {
    exchange_candidate_with_boundaries_and_external_heat(
        parcels,
        transfers,
        heat,
        boundaries,
        &[],
        endpoints,
        controls,
        cancelled,
    )
}

/// Joint internal transport, external phase streams, and additional heat, with
/// one PR endpoint reconstruction per affected node. The external owner must
/// accept its opposite mass/energy change atomically; no phase law or external
/// storage is synthesized. Boundary owner IDs and external-heat IDs are each
/// unique within their own explicitly separate physical channel.
#[allow(clippy::too_many_arguments)]
pub fn exchange_candidate_with_boundaries_and_external_heat(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
    transfers: &[SpeciesNetworkTransfer],
    heat: &[SpeciesNetworkHeatTransfer],
    boundaries: &[RealGasSpeciesNetworkBoundary<'_>],
    external_heat: &[SpeciesNetworkExternalHeat],
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<RealGasSpeciesNetworkCandidate, Error> {
    allowed_candidate(
        parcels,
        transfers,
        heat,
        boundaries,
        external_heat,
        &BTreeMap::new(),
        endpoints,
        controls,
        cancelled,
    )
}

// Every public entry lowers here; only the molecular-diffusion entry carries a
// nonempty mass allowance.
#[allow(clippy::too_many_arguments)]
fn allowed_candidate(
    parcels: &BTreeMap<u64, IsochoricRealGasParcel>,
    transfers: &[SpeciesNetworkTransfer],
    heat: &[SpeciesNetworkHeatTransfer],
    boundaries: &[RealGasSpeciesNetworkBoundary<'_>],
    external_heat: &[SpeciesNetworkExternalHeat],
    mass_allowance_kg: &BTreeMap<u64, f64>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RealGasSpeciesNetworkCandidate, Error> {
    check(&mut cancelled)?;
    snapshot::admit_work(
        parcels.len(),
        transfers.len(),
        heat.len(),
        external_heat.len(),
        boundary::work_count(boundaries)?,
    )?;
    let borrowed = parcels.iter().map(|(&id, parcel)| (id, parcel)).collect();
    exchange_borrowed_with_boundaries(
        &borrowed,
        transfers,
        heat,
        external_heat,
        boundaries,
        mass_allowance_kg,
        endpoints,
        controls,
        &mut cancelled,
    )
}

// The pair API lowers into this SAME owner without cloning its two complete
// material states merely to construct a temporary input map.
#[allow(clippy::too_many_arguments)]
pub(super) fn exchange_borrowed(
    parcels: &Parcels<'_>,
    transfers: &[SpeciesNetworkTransfer],
    heat: &[SpeciesNetworkHeatTransfer],
    external_heat: &[SpeciesNetworkExternalHeat],
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RealGasSpeciesNetworkCandidate, Error> {
    exchange_borrowed_with_boundaries(
        parcels,
        transfers,
        heat,
        external_heat,
        &[],
        &BTreeMap::new(),
        endpoints,
        controls,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn exchange_borrowed_with_boundaries(
    parcels: &Parcels<'_>,
    transfers: &[SpeciesNetworkTransfer],
    heat: &[SpeciesNetworkHeatTransfer],
    external_heat: &[SpeciesNetworkExternalHeat],
    boundaries: &[RealGasSpeciesNetworkBoundary<'_>],
    mass_allowance_kg: &BTreeMap<u64, f64>,
    endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: SpeciesNetworkSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RealGasSpeciesNetworkCandidate, Error> {
    let mut snapshot = snapshot::prepare(
        parcels,
        transfers,
        heat,
        external_heat,
        boundaries,
        mass_allowance_kg,
        endpoints,
        controls,
        cancelled,
    )?;
    let fluxes = caloric::streams(parcels, &mut snapshot, cancelled)?;
    caloric::finish(parcels, snapshot, fluxes, endpoints, controls, cancelled)
}

#[cfg(test)]
mod tests;
