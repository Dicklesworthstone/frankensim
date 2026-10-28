//! Pressure restrictions, porous leakage and selective transport on one snapshot.
//!
//! Sample nozzle, Darcy, membrane, Maxwell-Stefan, dusty-gas and effusion laws without changing gas.
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
/// Coupled molecular, Knudsen and viscous transport through one porous path.
pub mod dusty_gas;
/// Finite gas inventories linked by independently authored porous throats.
pub mod pore_domain;
/// Bidirectional free-molecular species and energy transfer through thin apertures.
pub mod effusion;
mod sample;
pub use sample::GasTransportRates;
pub(crate) use sample::sample_transfers;

use std::collections::BTreeMap;
use super::{ClosedGasParcel, GasEnergySolve, SpeciesNetworkError, SpeciesNetworkSolve,
    SpeciesNetworkStep, SpeciesNetworkTransfer, MAX_SNAPSHOT_TRANSFERS, inventory};
use super::super::diffusion::network::{self as molecular, DiffusionNetworkError,
    DiffusionNetworkLink, DiffusionNetworkLinkStep, DiffusionNetworkSolve};
use super::super::super::nozzle::network::snapshot::{self as pressure,
    SnapshotNozzleError, SnapshotNozzleLink, SnapshotNozzleLinkStep};
use porous::{PorousNetworkError, PorousNetworkLink, PorousNetworkLinkStep};
use membrane::{MembraneNetworkError, MembraneNetworkLink, MembraneNetworkLinkStep};
use dusty_gas::{DustyGasNetworkError, DustyGasNetworkLink, DustyGasNetworkLinkStep};
use effusion::{EffusionNetworkError, EffusionNetworkLink, EffusionNetworkLinkStep};

/// All authored mechanisms read the same original parcel inventory.
///
/// Empty slices close that mechanism. Paths may connect the same volumes, but
/// must describe distinct exposed areas or explicitly decomposed mechanisms;
/// adding a nozzle and Darcy law for the same aperture would double its flux.
/// A dusty-gas path already includes molecular, Knudsen and viscous transport;
/// do not add those mechanisms again for its same exposed pore area.
#[derive(Clone, Copy, Debug, Default)]
pub struct GasTransportPaths<'a> {
    pub nozzles: &'a [SnapshotNozzleLink],
    pub diffusion: &'a [DiffusionNetworkLink<'a>],
    pub porous: &'a [PorousNetworkLink<'a>],
    pub membranes: &'a [MembraneNetworkLink<'a>],
    pub dusty_gas: &'a [DustyGasNetworkLink<'a>],
    pub effusion: &'a [EffusionNetworkLink<'a>],
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
    DustyGas(DustyGasNetworkError),
    Effusion(EffusionNetworkError),
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

impl From<DustyGasNetworkError> for GasTransportError {
    fn from(error: DustyGasNetworkError) -> Self {
        match error { DustyGasNetworkError::Cancelled => Self::Cancelled, other => Self::DustyGas(other) }
    }
}

impl From<EffusionNetworkError> for GasTransportError {
    fn from(error: EffusionNetworkError) -> Self {
        match error { EffusionNetworkError::Cancelled => Self::Cancelled, other => Self::Effusion(other) }
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
    /// Coupled fluxes in the stationary-solid frame, without a zero-sum constraint.
    pub dusty_gas: Vec<DustyGasNetworkLinkStep>,
    /// Both original molecular beams and their kinetic energy convention.
    pub effusion: Vec<EffusionNetworkLinkStep>,
    /// Partial enthalpies and separately identified molecular energy corrections
    /// are applied together before one endpoint reconstruction.
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
    /// Dusty-gas cross-diffusion can drive an individual species uphill; imposing
    /// separate partial-pressure ordering would reject physically coupled fluxes.
    /// Effusive thermal transpiration also has no pressure-order gate. Both
    /// mechanisms instead share the inventory and total endpoint entropy gates.
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
            || self.dusty_gas.iter().any(|row| row.parcels.iter().any(|id| !parcels.contains_key(id)))
            || self.effusion.iter().any(|row| row.parcels.iter().any(|id| !parcels.contains_key(id)))
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

/// Advance all declared pressure, porous, diffusive, membrane and molecular paths transactionally.
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
    let sampled = sample_transfers(parcels, paths, duration_s, endpoints, controls, &mut cancelled)?;
    let proposal = super::exchange_candidate_with_heat(parcels, &sampled.transfers, &sampled.heat,
        endpoints, controls.exchange, &mut cancelled)?;
    sampled.rates.admit_endpoint_driving_forces(proposal.parcels(), &mut cancelled)?;
    let entropy_change_j_k = molecular::admit_entropy(parcels, proposal.parcels(),
        controls.entropy_tolerance_j_k, controls.exchange.mass_tolerance_kg, &mut cancelled)?;
    check(&mut cancelled)?;
    let (parcels, exchange) = proposal.into_parts();
    Ok(GasTransportCandidate { parcels, step:sampled.rates.into_closed_step(exchange, entropy_change_j_k) })
}

#[cfg(test)]
mod tests;
