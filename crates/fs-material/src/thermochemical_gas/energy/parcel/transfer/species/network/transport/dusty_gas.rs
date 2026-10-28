//! Coupled molecular, Knudsen and viscous pore transport on one gas snapshot.
//!
//! A dusty-gas face already includes the medium's Darcy contribution. Author it
//! as one physical connection; adding another Darcy or nozzle edge for the same
//! pore path would count the pressure-driven transport twice.
//!
//! This adapter only samples signed species rates. The shared selective-transfer
//! owner spends the original species inventories and carries their original
//! partial enthalpies once. Cross-diffusion can drive a species against its own
//! partial-pressure gradient, and unequal Knudsen mobilities permit a net molar
//! flux. Neither a per-species pressure-order gate nor a zero-molar-flux
//! constraint belongs on this path. The enclosing finite update admits gross
//! withdrawals, conserved mass/energy and the combined endpoint entropy.

use std::collections::BTreeMap;

use fs_matdb::QueryPoint;

use super::super::SpeciesNetworkTransfer;
use super::super::super::diffusion::{
    BinaryDiffusionMaterial,
    dusty_gas::{DustyGasFaceSolve, DustyGasFlux, DustyGasMedium},
};
use super::super::super::super::snapshot::{
    MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS,
};
use crate::thermochemical_gas::energy::GasEnergyError;
use crate::thermochemical_gas::energy::parcel::ClosedGasParcel;

/// A rigid porous connection whose one constitutive law includes viscous flow.
#[derive(Clone, Copy, Debug)]
pub struct DustyGasNetworkLink<'a> {
    /// Reporting identity; never selects a sequential inventory-update order.
    pub id: u64,
    /// Oriented left/right finite volumes.
    pub parcels: [u64; 2],
    pub medium: &'a DustyGasMedium,
    /// Exposed area [m²]. Zero closes the face without material-source queries.
    pub area_m2: f64,
    /// Source-bound binary molecular diffusivities for the joint species catalog.
    pub pairs: &'a [BinaryDiffusionMaterial<'a>],
    /// External gas-source coordinates retained by the face's material queries.
    pub gas_context: &'a QueryPoint,
    pub solve: DustyGasFaceSolve,
}

/// Original face result before the shared finite-inventory update.
#[derive(Clone, Debug, PartialEq)]
pub struct DustyGasNetworkLinkStep {
    pub id: u64,
    pub parcels: [u64; 2],
    pub flux: DustyGasFlux,
}

/// A refusal leaves every original gas parcel unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum DustyGasNetworkError {
    Invalid(&'static str),
    Face { id: u64, error: GasEnergyError },
    Cancelled,
}

impl core::fmt::Display for DustyGasNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material dusty-gas network: {self:?}")
    }
}
impl std::error::Error for DustyGasNetworkError {}

use DustyGasNetworkError as Error;

fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}

/// Sample every face from the unchanged snapshot, in canonical reporting order.
/// No accepted gas is mutated, and no species or enthalpy is applied here.
pub(super) fn sample_transfers(
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    links: &[DustyGasNetworkLink<'_>],
    duration_s: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(Vec<SpeciesNetworkTransfer>, Vec<DustyGasNetworkLinkStep>), Error> {
    check(&mut cancelled)?;
    if !duration_s.is_finite() || duration_s <= 0.0
        || parcels.is_empty() || parcels.len() > MAX_SNAPSHOT_PARCELS
        || links.len() > MAX_SNAPSHOT_TRANSFERS
    {
        return Err(Error::Invalid(
            "bounded nonempty dusty-gas graph and positive duration required",
        ));
    }
    let mut ordered = BTreeMap::new();
    for link in links {
        check(&mut cancelled)?;
        if link.parcels[0] == link.parcels[1]
            || link.parcels.iter().any(|id| !parcels.contains_key(id))
            || !link.area_m2.is_finite() || link.area_m2 < 0.0
            || ordered.insert(link.id, link).is_some()
        {
            return Err(Error::Invalid(
                "unique dusty-gas IDs, distinct existing endpoints and nonnegative area required",
            ));
        }
    }
    let mut transfers = Vec::new();
    let mut reports = Vec::with_capacity(links.len());
    for (id, link) in ordered {
        check(&mut cancelled)?;
        let [left, right] = link.parcels;
        let flux = parcels[&left].dusty_gas_flux_to(
            &parcels[&right], link.medium, link.area_m2, link.pairs,
            link.gas_context, link.solve, &mut cancelled,
        ).map_err(|error| if error == GasEnergyError::Cancelled {
            Error::Cancelled
        } else {
            Error::Face { id, error }
        })?;
        for species in &flux.species {
            check(&mut cancelled)?;
            let rate = species.mass_rate_left_to_right_kg_s;
            let mass = rate * duration_s;
            if !mass.is_finite() || (rate != 0.0 && mass == 0.0) {
                return Err(Error::Invalid(
                    "integrated dusty-gas species mass is unrepresentable",
                ));
            }
            if mass == 0.0 { continue; }
            if transfers.len() == MAX_SNAPSHOT_TRANSFERS {
                return Err(Error::Invalid(
                    "dusty-gas graph exceeds the selective-transfer bound",
                ));
            }
            let (donor, receiver) = if mass > 0.0 {
                (left, right)
            } else {
                (right, left)
            };
            transfers.push(SpeciesNetworkTransfer {
                donor, receiver, species: species.species.clone(), mass_kg: mass.abs(),
            });
        }
        reports.push(DustyGasNetworkLinkStep { id, parcels: link.parcels, flux });
    }
    check(&mut cancelled)?;
    Ok((transfers, reports))
}
