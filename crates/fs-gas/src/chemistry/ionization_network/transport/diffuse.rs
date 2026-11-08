//! Spectral diffuse radiation shared by co-located atomic charge-state ladders.
//!
//! Every source emits its original Milne spectrum; every admitted bound-free
//! channel competes for those photons through the same total opacity. Local
//! recycling between different species retains the difference in their chemical
//! thresholds as electron energy. No mean emitted photon, on-the-spot fraction,
//! excited-state cascade or additional atomic source is inferred.

use super::*;
use crate::chemistry::diffuse_recombination::evolution::DiffuseRecombinationExposure;
use crate::chemistry::diffuse_recombination::{
    DiffuseRecombinationBalance, DiffuseRecombinationCellBalance, DiffuseRecombinationError,
    DiffuseRecombinationFlux, DiffuseRecombinationSolve,
};
use fs_qty::Time;

pub(super) mod field;

/// Original network/beam controls and independent spectral/angular allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkDiffuseSolve {
    pub network: NetworkRadiationSolve,
    pub radiation: DiffuseRecombinationSolve,
}
impl NetworkDiffuseSolve {
    pub const fn new(network: NetworkRadiationSolve, radiation: DiffuseRecombinationSolve) -> Self {
        Self { network, radiation }
    }
}

/// Accepted spontaneous emission and recycling in one evolving mixture.
/// Counts are gross except the explicitly signed net moments. Same-source,
/// cross-source and nonlocal absorption are disjoint; boundary beams retain
/// their existing independent band accounts on the parent result.
/// The signed net moments include absorption from those boundary beams.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkDiffuseCellExposure {
    pub emitted_photons: f64,
    pub self_absorbed_photons: f64,
    pub cross_channel_absorbed_photons: f64,
    pub nonlocal_absorbed_photons: f64,
    pub diffuse_absorbed_energy_j: f64,
    pub emitted_energy_j: f64,
    pub captured_electron_energy_j: f64,
    pub net_ionization_extent_mol: f64,
    pub net_electron_energy_j: f64,
    pub net_gas_energy_j: f64,
}
impl NetworkDiffuseCellExposure {
    pub fn diffuse_absorbed_photons(&self) -> f64 {
        self.self_absorbed_photons
            + self.cross_channel_absorbed_photons
            + self.nonlocal_absorbed_photons
    }
}

/// Independent time-integrated boundary and gas-source conservation audits.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkDiffuseEvolutionBalance {
    pub emitted_photons: f64,
    pub emitted_energy_j: f64,
    pub diffuse_absorbed_photons: f64,
    pub diffuse_absorbed_energy_j: f64,
    pub incident_photons: f64,
    pub incident_energy_j: f64,
    pub escaped_photons: f64,
    pub escaped_energy_j: f64,
    /// Additional escaped energy from the explicitly thin free-free closure.
    /// It is not part of the bound-free directional interfaces or photon count.
    pub free_free_escaped_energy_j: f64,
    pub net_ionization_extent_mol: f64,
    pub net_radiation_energy_j: f64,
    /// N_A times the direct net radiative extent plus escape minus incidence.
    pub photon_residual: f64,
    /// Direct signed gas radiation source plus escape minus incidence.
    pub radiation_energy_residual_j: f64,
    /// Actual gas storage plus bound-free and free-free escape, minus incidence
    /// and external work.
    pub energy_residual_j: f64,
    /// Allowance confined to subtraction of bright incoming/outgoing energy.
    /// The independent gas-source audit keeps its original absolute tolerance.
    pub boundary_energy_roundoff_j: f64,
    pub radiation_evaluations: usize,
    pub radiation_work: usize,
    pub maximum_normalized_spectral_estimate: f64,
    pub maximum_normalized_angular_estimate: f64,
    pub maximum_normalized_tail_bound: f64,
}

/// Directional interface exposures and original controls for the live field.
#[derive(Clone, Debug, PartialEq)]
pub struct NetworkDiffuseEvolutionSummary {
    pub(super) interfaces: Vec<DiffuseRecombinationExposure>,
    pub(super) balance: NetworkDiffuseEvolutionBalance,
    pub(super) controls: NetworkDiffuseSolve,
}
impl NetworkDiffuseEvolutionSummary {
    pub fn interfaces(&self) -> &[DiffuseRecombinationExposure] {
        &self.interfaces
    }
    pub const fn balance(&self) -> NetworkDiffuseEvolutionBalance {
        self.balance
    }
    pub const fn controls(&self) -> NetworkDiffuseSolve {
        self.controls
    }
}

/// Evolve shared atomic ladders, temperatures and the diffuse radiation field.
///
/// The original collisional/thermal network surrounds a simultaneous midpoint
/// radiation update. Recombination photons may escape, recycle locally across
/// different thresholds, or ionize another cell. Every trial recomputes the
/// spectrum and opacity from its actual populations and temperature. Source
/// tables and transport buffers persist; no artificial electron seed, clipped
/// population, empirical escape fraction or mean photon replaces the sources.
///
/// The existing full/two-half controller also checks integrated radiation
/// moments and all interface exposures. Spectral/angular errors and omitted
/// tails have separate bounds, and all rejected work shares the global budget.
/// This quasi-steady plane-parallel closure does not model photon travel time,
/// scattering, excited-state cascades or unrestricted stiff kinetics.
pub fn advance_network_diffuse_recombination(
    cells: &[NetworkRadiationCell<'_>],
    area: Area,
    beams: &[DiffuseRecombinationBeam],
    duration: Time,
    controls: NetworkDiffuseSolve,
    cancelled: impl FnMut() -> bool,
) -> Result<super::NetworkRadiationResult, NetworkRadiationError> {
    super::evolution::advance_impl(
        cells,
        area,
        beams,
        duration,
        controls.network,
        Some(controls.radiation),
        cancelled,
    )
}

/// Spectral rates of one canonical source edge in its actual mixed cell.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkDiffuseChannelRate {
    pub emitted_photons_per_s: f64,
    pub emitted_power_w: f64,
    pub captured_electron_power_w: f64,
    /// Photons emitted and reabsorbed by this same original channel.
    pub self_absorbed_photons_per_s: f64,
    /// Photons emitted by another channel in the same mixed cell.
    pub cross_channel_absorbed_photons_per_s: f64,
    /// Photons emitted in other cells.
    pub nonlocal_absorbed_photons_per_s: f64,
    pub incident_absorbed_photons_per_s: f64,
    pub incident_absorbed_power_w: f64,
    /// Absorbed spontaneous radiation, retained separately from bright beams.
    pub diffuse_absorbed_power_w: f64,
    pub absorbed_power_w: f64,
    pub photoelectron_power_w: f64,
    pub net_ionization_rate_mol_per_s: f64,
    pub net_electron_power_w: f64,
    pub chemical_storage_rate_w: f64,
    pub net_gas_power_w: f64,
    pub emission_cutoff_j: f64,
}

/// Frozen source-retained mixture and independent cell/channel radiation rates.
#[derive(Clone, Debug, PartialEq)]
pub struct NetworkDiffuseFieldCell {
    initial: FrozenIonizedGasStateV1,
    sources: Vec<AtomicRadiativeRecombinationV1>,
    channels: Vec<NetworkDiffuseChannelRate>,
    beams: Vec<NetworkRadiationBeamRate>,
    balance: DiffuseRecombinationCellBalance,
}
impl NetworkDiffuseFieldCell {
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    pub fn sources(&self) -> &[AtomicRadiativeRecombinationV1] {
        &self.sources
    }
    pub fn channel_rates(&self) -> &[NetworkDiffuseChannelRate] {
        &self.channels
    }
    pub fn beam_rates(&self) -> &[NetworkRadiationBeamRate] {
        &self.beams
    }
    /// `self_absorbed_photons_per_s` here includes all same-cell recycling;
    /// channel rates distinguish same-source from cross-source absorption.
    /// Net rates are accumulated directly before per-channel f64 rounding.
    pub const fn balance(&self) -> DiffuseRecombinationCellBalance {
        self.balance
    }
}

/// Complete frozen field including source-support and quadrature diagnostics.
#[derive(Clone, Debug, PartialEq)]
pub struct NetworkDiffuseField {
    cells: Vec<NetworkDiffuseFieldCell>,
    beams: Vec<DiffuseRecombinationBeam>,
    beam_rates: Vec<NetworkRadiationBeamRate>,
    interfaces: Vec<DiffuseRecombinationFlux>,
    balance: DiffuseRecombinationBalance,
    area: Area,
    controls: NetworkDiffuseSolve,
}
impl NetworkDiffuseField {
    pub fn cells(&self) -> &[NetworkDiffuseFieldCell] {
        &self.cells
    }
    pub fn beams(&self) -> &[DiffuseRecombinationBeam] {
        &self.beams
    }
    pub fn beam_rates(&self) -> &[NetworkRadiationBeamRate] {
        &self.beam_rates
    }
    pub fn interfaces(&self) -> &[DiffuseRecombinationFlux] {
        &self.interfaces
    }
    pub const fn balance(&self) -> DiffuseRecombinationBalance {
        self.balance
    }
    pub const fn area(&self) -> Area {
        self.area
    }
    pub const fn controls(&self) -> NetworkDiffuseSolve {
        self.controls
    }
}

/// Evaluate original recombination spectra and shared-opacity recycling.
///
/// Local emission is isotropic and each mixed slab uses exact uniform-source
/// characteristics. All emitting and absorbing sources must support the common
/// spectrum; omitted tails must satisfy the original conditional bounds. In
/// particular, an absent or non-emitting channel cannot silently acquire an
/// unsupported high-energy opacity. This entry does not evolve the input state.
pub fn evaluate_network_diffuse_recombination(
    cells: &[NetworkRadiationCell<'_>],
    area: Area,
    beams: &[DiffuseRecombinationBeam],
    controls: NetworkDiffuseSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<NetworkDiffuseField, NetworkRadiationError> {
    controls.network.validate(cells, beams.len())?;
    poll(&mut cancelled)?;
    let problems = prepare_problems(cells, controls.network, &mut cancelled)?;
    let points: Vec<_> = problems.iter().map(Problem::initial_point).collect();
    let mut field =
        field::PreparedDiffuseField::new(cells, &problems, area, beams, controls, &mut cancelled)?;
    let remaining = controls
        .network
        .maximum_radiation_work
        .checked_sub(field.preparation_work())
        .ok_or(Error::WorkBudget)?;
    field.evaluate(&points, remaining, &mut cancelled)?;
    let mut results = reserved(cells.len())?;
    for (c, cell) in cells.iter().enumerate() {
        poll(&mut cancelled)?;
        results.push(NetworkDiffuseFieldCell {
            initial: cell.initial.clone(),
            sources: problems[c].sources().cloned().collect(),
            channels: field.channel_rates(c).to_vec(),
            beams: field.cell_beam_rates(c).to_vec(),
            balance: field.cell_balance(c),
        });
    }
    poll(&mut cancelled)?;
    Ok(NetworkDiffuseField {
        cells: results,
        beams: beams.to_vec(),
        beam_rates: field.beam_rates.clone(),
        interfaces: field.interfaces.clone(),
        balance: field.balance,
        area,
        controls,
    })
}

pub(super) fn diffuse_error(error: DiffuseRecombinationError) -> NetworkRadiationError {
    if matches!(error, DiffuseRecombinationError::Cancelled) {
        NetworkRadiationError::Cancelled
    } else {
        NetworkRadiationError::Diffuse(error)
    }
}
