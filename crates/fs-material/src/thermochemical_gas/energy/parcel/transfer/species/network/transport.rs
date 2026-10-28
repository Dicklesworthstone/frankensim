//! Pressure restrictions, porous leakage and selective transport on one snapshot.
//!
//! Sample nozzle, Darcy, membrane and Maxwell-Stefan laws without changing gas.
//! Lower bulk outflow to its original species masses, combine it with selective
//! fluxes, then invoke the ONE conservative species/enthalpy update. No inlet can
//! replenish an outgoing allowance or forward arrivals in this explicit step.
//!
//! Connection lists declare separate physical paths, not two models of
//! the same aperture. In particular, a molar-average diffusion velocity is not
//! silently equated to a nozzle's bulk velocity. This is not a spatial CFD
//! advection-diffusion discretization, pressure diffusion, or jet momentum law.

/// Isothermal pressure-driven leakage through actual porous materials.
pub mod porous;
/// Species-selective solution-diffusion through dense membranes.
pub mod membrane;

use std::collections::BTreeMap;
use super::{ClosedGasParcel, GasEnergySolve, SpeciesNetworkError, SpeciesNetworkSolve,
    SpeciesNetworkStep, SpeciesNetworkTransfer, MAX_SNAPSHOT_TRANSFERS, inventory};
use super::super::diffusion::network::{self as molecular, DiffusionNetworkError,
    DiffusionNetworkLink, DiffusionNetworkLinkStep, DiffusionNetworkSolve};
use super::super::super::nozzle::network::snapshot::{self as pressure,
    SnapshotNozzleError, SnapshotNozzleLink, SnapshotNozzleLinkStep};
use porous::{PorousNetworkError, PorousNetworkLink, PorousNetworkLinkStep};
use membrane::{MembraneNetworkError, MembraneNetworkLink, MembraneNetworkLinkStep};

/// All authored mechanisms read the same original parcel inventory.
///
/// Empty slices close that mechanism. Paths may connect the same volumes, but
/// must describe distinct exposed areas or explicitly decomposed mechanisms;
/// adding a nozzle and Darcy law for the same aperture would double its flux.
#[derive(Clone, Copy, Debug, Default)]
pub struct GasTransportPaths<'a> {
    pub nozzles: &'a [SnapshotNozzleLink],
    pub diffusion: &'a [DiffusionNetworkLink<'a>],
    pub porous: &'a [PorousNetworkLink<'a>],
    pub membranes: &'a [MembraneNetworkLink<'a>],
}

/// One aggregate allowance for all original species, irrespective of mechanism.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasTransportSolve {
    /// Original species owner checks gross withdrawal, local/global mass and energy.
    pub exchange: SpeciesNetworkSolve,
    /// Allowed rounding decrease in total closed-network endpoint entropy [J/K].
    pub entropy_tolerance_j_k: f64,
}

/// Sampling, reconstruction or combined endpoint refusal; inputs are immutable.
#[derive(Clone, Debug, PartialEq)]
pub enum GasTransportError {
    Invalid(&'static str),
    Nozzle(SnapshotNozzleError),
    Diffusion(DiffusionNetworkError),
    Porous(PorousNetworkError),
    Membrane(MembraneNetworkError),
    Exchange(SpeciesNetworkError),
    Cancelled,
}
impl core::fmt::Display for GasTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "combined material gas transport: {self:?}")
    }
}
impl std::error::Error for GasTransportError {}
impl From<SnapshotNozzleError> for GasTransportError {
    fn from(error: SnapshotNozzleError) -> Self {
        match error { SnapshotNozzleError::Cancelled => Self::Cancelled, other => Self::Nozzle(other) }
    }
}
impl From<DiffusionNetworkError> for GasTransportError {
    fn from(error: DiffusionNetworkError) -> Self {
        match error { DiffusionNetworkError::Cancelled => Self::Cancelled, other => Self::Diffusion(other) }
    }
}
impl From<SpeciesNetworkError> for GasTransportError {
    fn from(error: SpeciesNetworkError) -> Self {
        match error { SpeciesNetworkError::Cancelled => Self::Cancelled, other => Self::Exchange(other) }
    }
}
impl From<PorousNetworkError> for GasTransportError {
    fn from(error: PorousNetworkError) -> Self {
        match error { PorousNetworkError::Cancelled => Self::Cancelled, other => Self::Porous(other) }
    }
}
impl From<MembraneNetworkError> for GasTransportError {
    fn from(error: MembraneNetworkError) -> Self {
        match error { MembraneNetworkError::Cancelled => Self::Cancelled, other => Self::Membrane(other) }
    }
}

/// Original source-bearing rates plus the single accepted inventory accounting.
#[derive(Clone, Debug, PartialEq)]
pub struct GasTransportStep {
    pub duration_s: f64,
    /// Nozzle IDs are unique within this list; they do not order updates.
    pub nozzles: Vec<SnapshotNozzleLinkStep>,
    /// Diffusion IDs are unique within this list; they do not order updates.
    pub diffusion: Vec<DiffusionNetworkLinkStep>,
    /// Source-resolved porous rates before any accepted inventory was changed.
    pub porous: Vec<PorousNetworkLinkStep>,
    /// Per-species partial-pressure-driven permeation on the same snapshot.
    pub membranes: Vec<MembraneNetworkLinkStep>,
    /// Partial enthalpies, including formation references, are transported once.
    pub exchange: SpeciesNetworkStep,
    /// Conserved-species entropy change; species-mass residuals stay separate.
    /// This is not a stability/accuracy certificate.
    pub entropy_change_j_k: f64,
}

impl GasTransportStep {
    /// Reapply ORIGINAL transport driving-force ordering after an enclosing
    /// thermal/mechanical proposal changes the transported gas endpoints.
    ///
    /// This does not sample another flow or spend inventory again. The supplied
    /// map must retain every original parcel identity; the actual positive
    /// transfers and directions in this report remain authoritative. Nozzles and
    /// porous paths retain total-pressure order; membranes retain each active
    /// species' partial-pressure order. Heat/work can change entropy, and the transport
    /// entropy report is not a certificate for that subsequent process.
    pub fn admit_endpoint_nozzle_pressures(
        &self, parcels: &BTreeMap<u64, ClosedGasParcel>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), GasTransportError> {
        self.admit_endpoint_driving_forces(parcels, cancelled)
    }

    /// Admit every sampled pressure/permeation mechanism on the actual endpoint.
    /// The original nozzle-named entry point delegates here so existing coupled
    /// consumers cannot accidentally bypass newly authored porous or membrane paths.
    pub fn admit_endpoint_driving_forces(
        &self, parcels: &BTreeMap<u64, ClosedGasParcel>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), GasTransportError> {
        check(&mut cancelled)?;
        if !parcels.keys().copied().eq(self.exchange.parcels.iter().map(|row| row.id))
            || self.nozzles.iter().any(|row| !parcels.contains_key(&row.donor)
                || !parcels.contains_key(&row.receiver))
            || self.porous.iter().any(|row| row.parcels.iter().any(|id| !parcels.contains_key(id)))
            || self.membranes.iter().any(|row| row.parcels.iter().any(|id| !parcels.contains_key(id)))
        {
            return Err(GasTransportError::Invalid(
                "final transport admission needs the complete original parcel identity map",
            ));
        }
        pressure::admit_pressures(parcels, &self.nozzles, &mut cancelled)?;
        porous::admit_pressures(parcels, &self.porous, &mut cancelled)?;
        membrane::admit_partial_pressures(parcels, &self.membranes, &mut cancelled)?;
        check(&mut cancelled)
    }
}

/// All proposed material states and their common balances, or no endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct GasTransportCandidate {
    parcels: BTreeMap<u64, ClosedGasParcel>,
    step: GasTransportStep,
}
impl GasTransportCandidate {
    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> { &self.parcels }
    pub const fn step(&self) -> &GasTransportStep { &self.step }
    pub fn into_parts(self) -> (BTreeMap<u64, ClosedGasParcel>, GasTransportStep) {
        (self.parcels, self.step)
    }
}
use GasTransportError as Error;
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}

/// Advance explicitly declared pressure-driven AND diffusive connections once.
///
/// All laws sample the same unchanged volumes, compositions, temperatures and
/// pressures. Existing samplers retain their source-domain, throat and face
/// admission. The selective owner enforces one GROSS withdrawal allowance per
/// original species across BOTH mechanisms before one reconstruction per volume.
/// Thus individually safe paths may be refused together; neither incoming flux
/// nor an opposite transfer can conceal excessive gross withdrawal.
///
/// The combined endpoints must retain every initially active nozzle's pressure
/// ordering and satisfy the existing extensive entropy gate. Diffusion remains
/// near-isothermal/isobaric in the molar-average frame on its own declared paths;
/// adjacent nozzle pressure drops do not relax those assumptions. Parallel paths
/// between the same volumes are allowed only as distinct authored mechanisms.
/// There is no extra flow work, chemical heat, inferred diffusivity or volume
/// update. Caller-controlled timestep refinement remains required (first order).
///
/// Empty lists close that mechanism. Both empty preserve all gases exactly.
/// No previously received species can leave another connection in this step.
pub fn step_candidate(
    parcels: &BTreeMap<u64, ClosedGasParcel>, nozzles: &[SnapshotNozzleLink],
    diffusion: &[DiffusionNetworkLink<'_>], duration_s: f64,
    endpoints: &BTreeMap<u64, GasEnergySolve>, controls: GasTransportSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<GasTransportCandidate, Error> {
    step_candidate_with_paths(parcels, GasTransportPaths {
        nozzles, diffusion, ..GasTransportPaths::default()
    }, duration_s, endpoints, controls, cancelled)
}

/// Advance all pressure, porous, diffusive and membrane paths transactionally.
///
/// Every mechanism samples the original gases. Their combined GROSS species
/// withdrawal is admitted before one caloric reconstruction per finite volume.
/// Pressure-order and entropy checks apply to the FINAL combined endpoint.
/// The caller controls first-order time refinement; failed proposals publish
/// no partial state. Constitutive source domains remain owned by each law.
pub fn step_candidate_with_paths(
    parcels: &BTreeMap<u64, ClosedGasParcel>, paths: GasTransportPaths<'_>,
    duration_s: f64, endpoints: &BTreeMap<u64, GasEnergySolve>,
    controls: GasTransportSolve, mut cancelled: impl FnMut() -> bool,
) -> Result<GasTransportCandidate, Error> {
    check(&mut cancelled)?;
    let count = [paths.nozzles.len(), paths.diffusion.len(), paths.porous.len(), paths.membranes.len()]
        .into_iter().try_fold(0usize, |total, n| total.checked_add(n));
    if count.is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS) {
        return Err(Error::Invalid("combined connection count exceeds the shared transport bound"));
    }
    let (bulk, nozzle_reports) = pressure::sample_transfers(
        parcels, paths.nozzles, duration_s, endpoints, &mut cancelled)?;
    let (mut transfers, diffusion_reports) = molecular::sample_transfers(
        parcels, paths.diffusion, duration_s, endpoints,
        DiffusionNetworkSolve { exchange: controls.exchange,
            entropy_tolerance_j_k: controls.entropy_tolerance_j_k }, &mut cancelled)?;
    let (porous_bulk, porous_reports) = porous::sample_transfers(
        parcels, paths.porous, duration_s, &mut cancelled)?;
    let (selective, membrane_reports) = membrane::sample_transfers(
        parcels, paths.membranes, duration_s, &mut cancelled)?;
    if transfers.len().checked_add(selective.len()).is_none_or(|n| n > MAX_SNAPSHOT_TRANSFERS) {
        return Err(Error::Invalid("combined selective flux count exceeds the shared work bound"));
    }
    transfers.extend(selective);
    for row in bulk.into_iter().chain(porous_bulk) {
        check(&mut cancelled)?;
        if row.mass_kg == 0.0 { continue; }
        // The original composition owner converts mole fractions to mass.
        // Summing dm_i*h_i is exactly the bulk upstream dm*h mixture law;
        // the shared species owner evaluates h_i with original source references.
        let parts = inventory(parcels[&row.donor].model(), row.mass_kg)
            .map_err(|error| super::material(row.donor, error))?;
        for (species, mass_kg) in parts {
            check(&mut cancelled)?;
            if transfers.len() == MAX_SNAPSHOT_TRANSFERS {
                return Err(Error::Invalid("combined species flux count exceeds the shared work bound"));
            }
            transfers.push(SpeciesNetworkTransfer { donor: row.donor,
                receiver: row.receiver, species, mass_kg });
        }
    }
    let proposal = super::exchange_candidate(parcels, &transfers, endpoints,
        controls.exchange, &mut cancelled)?;
    // Evaluate pressure crossing on the FINAL combined inventory, not on either
    // mechanism in isolation. Neither failed proposal publishes a subset.
    pressure::admit_pressures(proposal.parcels(), &nozzle_reports, &mut cancelled)?;
    porous::admit_pressures(proposal.parcels(), &porous_reports, &mut cancelled)?;
    membrane::admit_partial_pressures(proposal.parcels(), &membrane_reports, &mut cancelled)?;
    let entropy_change_j_k = molecular::admit_entropy(parcels, proposal.parcels(),
        controls.entropy_tolerance_j_k, controls.exchange.mass_tolerance_kg, &mut cancelled)?;
    check(&mut cancelled)?;
    let (parcels, exchange) = proposal.into_parts();
    Ok(GasTransportCandidate { parcels, step: GasTransportStep { duration_s,
        nozzles: nozzle_reports, diffusion: diffusion_reports, porous: porous_reports,
        membranes: membrane_reports, exchange, entropy_change_j_k } })
}

#[cfg(test)]
mod tests;

