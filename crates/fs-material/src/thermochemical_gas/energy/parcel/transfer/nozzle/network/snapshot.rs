//! Unsplit pressure-driven exchange: every nozzle sees the SAME accepted state.
//!
//! Resolve original NASA-9 throat flows first, then delegate all integrated
//! species/enthalpy fluxes to the shared multiport material update. One new
//! thermodynamic state is reconstructed per affected volume, not per link.
//! Link ids affect only reporting, never flow arithmetic or source histories.
//!
//! This is an explicit first-order graph step, not an implicit network solve.
//! Original quasi-steady/choked-flow and well-mixed fixed-volume limits remain.
//! Closed/equilibrated links are zero at the initial snapshot; new pressure
//! differences drive them on the next interval. Refine time at events.

use std::collections::BTreeMap;
use super::super::{GasNozzle, GasNozzleFlow, GasEnergySolve, GasEnergyError, ClosedGasParcel};
use super::super::super::snapshot::{self, GasSnapshotError, GasSnapshotSolve,
    GasSnapshotStep, GasSnapshotTransfer, MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS};

/// Physical restriction and its own throat-source admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SnapshotNozzleLink {
    /// Reporting identity; does not choose the update order.
    pub id: u64,
    /// Two distinct existing volumes. Upstream follows their initial pressure.
    pub parcels: [u64; 2],
    /// Original area and discharge coefficient; zero closes the connection.
    pub nozzle: GasNozzle,
    /// Original caloric bracket, tolerance and work bound for throat resolution.
    pub nozzle_solve: GasEnergySolve,
}
/// Actual original-snapshot flow retained with its source-bearing throat state.
#[derive(Clone, Debug, PartialEq)]
pub struct SnapshotNozzleLinkStep {
    /// Original link identity.
    pub id: u64,
    /// Original upstream volume; at equal pressure the smaller id is selected.
    pub donor: u64,
    /// Original receiving volume.
    pub receiver: u64,
    /// Existing nozzle owner's complete result.
    pub flow: GasNozzleFlow,
    /// Integrated original rate over the physical duration [kg].
    pub transferred_mass_kg: f64,
}
/// Every physical account belongs to this one all-volume interval.
#[derive(Clone, Debug, PartialEq)]
pub struct SnapshotNozzleStep {
    /// Actual physical duration [s].
    pub duration_s: f64,
    /// Source-resolved flows in reporting-id order.
    pub links: Vec<SnapshotNozzleLinkStep>,
    /// Actual shared-state local/global species, mass and energy accounts.
    pub exchange: GasSnapshotStep,
}
/// Full graph proposal. Accepted input parcels never change during iteration.
#[derive(Clone, Debug, PartialEq)]
pub struct SnapshotNozzleCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: SnapshotNozzleStep,
}
impl SnapshotNozzleCandidate {
    /// All complete endpoints, ready for existing thermal/mechanical consumers.
    #[must_use]
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> { &self.parcels }
    /// Original material flows and checked endpoint balances.
    #[must_use]
    pub const fn step(&self) -> &SnapshotNozzleStep { &self.step }
    /// Publish every endpoint together, with no second heat or mass deposition.
    #[must_use]
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, SnapshotNozzleStep) {
        (self.parcels, self.step)
    }
}
/// A failure never returns some accepted graph nodes and some old nodes.
#[derive(Clone, Debug, PartialEq)]
pub enum SnapshotNozzleError {
    /// Invalid graph geometry, clock or numerical admission.
    Invalid(&'static str),
    /// The original throat owner refused this physical connection.
    Link { id: u64, error: GasEnergyError },
    /// Shared conservative multiport reconstruction refused.
    Transfer(GasSnapshotError),
    /// The combined graph step crossed an initially active pressure ordering.
    /// Shorten the physical interval; do not clip a mass to make it fit.
    EquilibriumCrossed { id: u64 },
    /// Cancelled before returning all endpoints.
    Cancelled,
}
impl core::fmt::Display for SnapshotNozzleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "snapshot nozzle network: {self:?}")
    }
}
impl std::error::Error for SnapshotNozzleError {}
use SnapshotNozzleError as Error;
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}

/// Evaluate every restriction against the same immutable physical network.
///
/// Unlike the existing Lie-split API, changing link IDs cannot change the
/// resulting state. All incoming/outgoing source inventories are accumulated
/// before any endpoint is reconstructed. `endpoints` assigns exactly one
/// original energy solver to each volume, independent of its changing role as
/// an inlet or outlet. Global controls also bound local balances and aggregate
/// GROSS withdrawal from the original inventory, not replenished inventory.
///
/// After joint reconstruction, every initially active connection must retain
/// its pressure ordering. Other connections can cause an equilibrium crossing
/// even if the same edge would have passed an isolated pair step. This check
/// refuses the entire interval. It is not an event locator or a time-error bound.
pub fn step_candidate(
    parcels: &BTreeMap<u64, ClosedGasParcel>, links: &[SnapshotNozzleLink], duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>, controls: GasSnapshotSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SnapshotNozzleCandidate, Error> {
    let (transfers, reports) = sample_transfers(parcels, links, duration_s, endpoints, &mut cancelled)?;
    let candidate = snapshot::transfer_candidate(parcels, &transfers, endpoints, controls, &mut cancelled)
        .map_err(|error| if error == GasSnapshotError::Cancelled { Error::Cancelled } else { Error::Transfer(error) })?;
    admit_pressures(candidate.parcels(), &reports, &mut cancelled)?;
    check(&mut cancelled)?;
    let (parcels, exchange) = candidate.into_parts();
    Ok(SnapshotNozzleCandidate { parcels, step: SnapshotNozzleStep { duration_s, links: reports, exchange } })
}

// Shared by the nozzle-only and combined transport paths. This evaluates rates
// on the original snapshot but performs NO inventory or endpoint update.
pub(crate) fn sample_transfers(
    parcels: &BTreeMap<u64, ClosedGasParcel>, links: &[SnapshotNozzleLink], duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>, mut cancelled: impl FnMut() -> bool,
) -> Result<(Vec<GasSnapshotTransfer>, Vec<SnapshotNozzleLinkStep>), Error> {
    check(&mut cancelled)?;
    if !duration_s.is_finite() || duration_s <= 0.0 || parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS || links.len() > MAX_SNAPSHOT_TRANSFERS
        || endpoints.len() != parcels.len() || !endpoints.keys().eq(parcels.keys()) {
        return Err(Error::Invalid("positive physical duration, bounded graph and exact endpoint solver map required"));
    }
    let mut ordered = BTreeMap::new();
    for link in links {
        check(&mut cancelled)?;
        if link.parcels[0] == link.parcels[1] || link.parcels.iter().any(|id| !parcels.contains_key(id))
            || ordered.insert(link.id, link).is_some() {
            return Err(Error::Invalid("unique reporting IDs and distinct existing link endpoints required"));
        }
    }
    let mut transfers = Vec::with_capacity(links.len());
    let mut reports = Vec::with_capacity(links.len());
    for (&id, link) in &ordered {
        check(&mut cancelled)?;
        let [a, b] = link.parcels;
        let (pa, pb) = (parcels[&a].state().state().pressure, parcels[&b].state().state().pressure);
        let (donor, receiver) = if pa > pb || (pa == pb && a < b) { (a, b) } else { (b, a) };
        let flow = parcels[&donor].nozzle_flow_to(parcels[&receiver].state().state().pressure,
            link.nozzle, link.nozzle_solve, &mut cancelled).map_err(|error| {
                if error == GasEnergyError::Cancelled { Error::Cancelled } else { Error::Link { id, error } }
            })?;
        let mass = flow.mass_flow_kg_s * duration_s;
        if !mass.is_finite() || (flow.mass_flow_kg_s > 0.0 && mass == 0.0) {
            return Err(Error::Invalid("integrated nozzle mass is unrepresentable"));
        }
        transfers.push(GasSnapshotTransfer { donor, receiver, mass_kg: mass });
        reports.push(SnapshotNozzleLinkStep { id, donor, receiver, flow, transferred_mass_kg: mass });
    }
    Ok((transfers, reports))
}

// The same pressure-order admission is applied to the JOINT transport endpoint.
pub(crate) fn admit_pressures(
    parcels: &BTreeMap<u64, ClosedGasParcel>, reports: &[SnapshotNozzleLinkStep],
    mut cancelled: impl FnMut() -> bool,
) -> Result<(), Error> {
    for report in reports {
        check(&mut cancelled)?;
        if report.transferred_mass_kg > 0.0
            && parcels[&report.donor].state().state().pressure
                < parcels[&report.receiver].state().state().pressure {
            return Err(Error::EquilibriumCrossed { id: report.id });
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
