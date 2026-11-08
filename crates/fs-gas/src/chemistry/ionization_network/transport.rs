//! Photon-conserving transport through co-located atomic charge-state ladders.
//!
//! Each cell has one original population vector and electron pool. At a photon
//! energy all its admitted absorbers compete through their summed optical
//! depth; absorption is partitioned in proportion to each channel's opacity.
//! Adjacent charge stages are not separate spatial layers. Source ordering is
//! canonical, and every reachable absorber's source support is checked even
//! when its initial population is zero. Cross sections prepare once.
//!
//! Slabs are uniformly mixed, fixed-volume and transversely infinite. Boundary
//! beams are discrete and quasi-static. Their number flux is per area normal
//! to the ray, so face throughput includes abs(mu) and optical depth includes
//! its reciprocal. This owner does not infer diffuse emission, scattering,
//! excited levels, secondary nonthermal ionization or photon flight time.

use super::solve::{Point, Problem};
use super::{IonizationNetworkError, IonizationNetworkSolve};
use crate::chemistry::collisional_radiative::CollisionalRadiativeThermal;
use crate::chemistry::diffuse_recombination::DiffuseRecombinationBeam;
use fs_qty::Area;
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::atomic_ionization::photoionization::{
    MAX_PHOTON_GROUPS_V1, PhotoionizationErrorV1,
};
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1;

mod field;

/// One original gas mixture and all explicitly modeled adjacent-stage sources.
#[derive(Clone, Copy, Debug)]
pub struct NetworkRadiationCell<'a> {
    pub initial: &'a FrozenIonizedGasStateV1,
    pub sources: &'a [AtomicRadiativeRecombinationV1],
    pub thermal: Option<&'a CollisionalRadiativeThermal>,
}
impl<'a> NetworkRadiationCell<'a> {
    pub const fn new(
        initial: &'a FrozenIonizedGasStateV1,
        sources: &'a [AtomicRadiativeRecombinationV1],
    ) -> Self {
        Self {
            initial,
            sources,
            thermal: None,
        }
    }
    pub const fn with_thermal(mut self, thermal: &'a CollisionalRadiativeThermal) -> Self {
        self.thermal = Some(thermal);
        self
    }
}

/// Source, temporal, allocation and total radiation-work bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkRadiationSolve {
    pub network: IonizationNetworkSolve,
    /// Roundoff allowance relative to incoming photons at each interface audit.
    pub relative_photon_tolerance: f64,
    pub maximum_cells: usize,
    /// Sum of all cell channel counts, including initially absent stages.
    pub maximum_cell_channels: usize,
    /// Beam times total-channel product, bounding the immutable source table.
    pub maximum_beam_channels: usize,
    /// Cell times beam product, bounding interface-rate and exposure buffers.
    pub maximum_cell_beams: usize,
    /// Includes predictor, accepted and rejected radiation queries.
    pub maximum_radiation_evaluations: usize,
    /// Includes source-table admission and all attempted characteristic sweeps.
    pub maximum_radiation_work: usize,
}
impl NetworkRadiationSolve {
    pub fn new(network: IonizationNetworkSolve) -> Self {
        Self {
            network,
            relative_photon_tolerance: 1e-12,
            maximum_cells: 4096,
            maximum_cell_channels: 65_536,
            maximum_beam_channels: 1_048_576,
            maximum_cell_beams: 1_048_576,
            maximum_radiation_evaluations: 1_000_000,
            maximum_radiation_work: 100_000_000,
        }
    }
    fn validate(self, cells: &[NetworkRadiationCell<'_>], beams: usize) -> Result<(), Error> {
        self.network.validate().map_err(Error::Chemistry)?;
        if !positive(self.relative_photon_tolerance)
            || self.relative_photon_tolerance >= 1.0
            || self.maximum_cells == 0
            || self.maximum_cells > 1_048_576
            || self.maximum_cell_channels == 0
            || self.maximum_cell_channels > 1_048_576
            || self.maximum_beam_channels == 0
            || self.maximum_beam_channels > 16_777_216
            || self.maximum_cell_beams == 0
            || self.maximum_cell_beams > 16_777_216
            || self.maximum_radiation_evaluations == 0
            || self.maximum_radiation_work == 0
        {
            return Err(Error::Invalid("network radiation controls"));
        }
        let channels = cells.iter().try_fold(0_usize, |n, cell| {
            if cell.sources.len() > self.network.maximum_channels {
                return Err(Error::WorkBudget);
            }
            n.checked_add(cell.sources.len()).ok_or(Error::WorkBudget)
        })?;
        if cells.len() > self.maximum_cells
            || beams > MAX_PHOTON_GROUPS_V1
            || channels > self.maximum_cell_channels
            || channels
                .checked_mul(beams)
                .is_none_or(|n| n > self.maximum_beam_channels)
            || cells
                .len()
                .checked_mul(beams)
                .is_none_or(|n| n > self.maximum_cell_beams)
        {
            return Err(Error::WorkBudget);
        }
        Ok(())
    }
}

/// Rate for one canonical source edge, after competition for shared photons.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkPhotoionizationRate {
    pub ionization_rate_mol_per_s: f64,
    pub absorbed_power_w: f64,
    pub electron_heating_power_w: f64,
}

/// One beam's local or whole-column photon and absorbed-energy rates.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetworkRadiationBeamRate {
    pub incident_photons_per_s: f64,
    pub absorbed_photons_per_s: f64,
    pub transmitted_photons_per_s: f64,
    pub absorbed_power_w: f64,
    /// Absorbed plus transmitted minus incident, accumulated before rounding.
    pub photon_residual_per_s: f64,
}

/// Frozen original mixture with rates in canonical source and input beam order.
#[derive(Clone, Debug, PartialEq)]
pub struct NetworkRadiationFieldCell {
    initial: FrozenIonizedGasStateV1,
    sources: Vec<AtomicRadiativeRecombinationV1>,
    channel_rates: Vec<NetworkPhotoionizationRate>,
    beam_rates: Vec<NetworkRadiationBeamRate>,
}
impl NetworkRadiationFieldCell {
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    pub fn sources(&self) -> &[AtomicRadiativeRecombinationV1] {
        &self.sources
    }
    pub fn channel_rates(&self) -> &[NetworkPhotoionizationRate] {
        &self.channel_rates
    }
    pub fn beam_rates(&self) -> &[NetworkRadiationBeamRate] {
        &self.beam_rates
    }
}

/// Source-retained frozen query. No chemistry or thermal evolution is implied.
#[derive(Clone, Debug, PartialEq)]
pub struct NetworkRadiationField {
    cells: Vec<NetworkRadiationFieldCell>,
    beams: Vec<DiffuseRecombinationBeam>,
    beam_rates: Vec<NetworkRadiationBeamRate>,
    area: Area,
    controls: NetworkRadiationSolve,
    radiation_work: usize,
}
impl NetworkRadiationField {
    pub fn cells(&self) -> &[NetworkRadiationFieldCell] {
        &self.cells
    }
    pub fn beams(&self) -> &[DiffuseRecombinationBeam] {
        &self.beams
    }
    pub fn beam_rates(&self) -> &[NetworkRadiationBeamRate] {
        &self.beam_rates
    }
    pub const fn area(&self) -> Area {
        self.area
    }
    pub const fn controls(&self) -> NetworkRadiationSolve {
        self.controls
    }
    pub const fn radiation_work(&self) -> usize {
        self.radiation_work
    }
}

/// Admission, support, conservation or work refusal without caller mutation.
#[derive(Clone, Debug, PartialEq)]
pub enum NetworkRadiationError {
    Invalid(&'static str),
    Chemistry(IonizationNetworkError),
    Cell {
        index: usize,
        error: IonizationNetworkError,
    },
    Source {
        cell: usize,
        channel: usize,
        beam: usize,
        error: PhotoionizationErrorV1,
    },
    WorkBudget,
    StepBudget,
    Arithmetic(&'static str),
    PhotonBalance {
        cell: Option<usize>,
        beam: usize,
        residual_photons: f64,
    },
    EnergyBalance {
        residual_j: f64,
    },
    Cancelled,
}
impl core::fmt::Display for NetworkRadiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "network photon transport: {self:?}")
    }
}
impl std::error::Error for NetworkRadiationError {}
type Error = NetworkRadiationError;

/// Evaluate shared-opacity beam attenuation at the original frozen populations.
///
/// At each energy, tau=sum(tau_channel), and each source absorbs its fraction
/// tau_channel/tau of the same absorbed photon rate. This preserves invariance
/// under relabelling identical absorbers (Pawlik & Schaye, arXiv:1008.1071,
/// Sec.4.1). A zero-opacity beam traverses the column unchanged. Every original
/// source must support every supplied beam, including a dark or absent channel.
pub fn evaluate_network_radiation_transport(
    cells: &[NetworkRadiationCell<'_>],
    area: Area,
    beams: &[DiffuseRecombinationBeam],
    controls: NetworkRadiationSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<NetworkRadiationField, Error> {
    controls.validate(cells, beams.len())?;
    let problems = prepare_problems(cells, controls, &mut cancelled)?;
    let points: Vec<_> = problems.iter().map(Problem::initial_point).collect();
    let mut field =
        field::PreparedField::new(cells, &problems, area, beams, controls, &mut cancelled)?;
    field.evaluate(&points, &mut cancelled)?;
    let mut results = reserved(cells.len())?;
    for (c, cell) in cells.iter().enumerate() {
        poll(&mut cancelled)?;
        results.push(NetworkRadiationFieldCell {
            initial: cell.initial.clone(),
            sources: problems[c].sources().cloned().collect(),
            channel_rates: field.channel_rates(c).to_vec(),
            beam_rates: field.cell_beam_rates(c).to_vec(),
        });
    }
    poll(&mut cancelled)?;
    Ok(NetworkRadiationField {
        cells: results,
        beams: beams.to_vec(),
        beam_rates: field.beam_rates.clone(),
        area,
        controls,
        radiation_work: field.work,
    })
}

fn prepare_problems<'a>(
    cells: &[NetworkRadiationCell<'a>],
    controls: NetworkRadiationSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Problem<'a>>, Error> {
    let mut problems = reserved(cells.len())?;
    for (c, cell) in cells.iter().enumerate() {
        poll(cancelled)?;
        problems.push(
            Problem::new(
                cell.initial,
                cell.sources,
                &[],
                cell.thermal,
                controls.network,
                cancelled,
            )
            .map_err(|error| cell_error(c, error))?,
        );
    }
    Ok(problems)
}
fn reserved<T>(n: usize) -> Result<Vec<T>, Error> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| Error::Arithmetic("bounded allocation"))?;
    Ok(values)
}
fn buffer<T: Clone>(n: usize, value: T) -> Result<Vec<T>, Error> {
    let mut values = reserved(n)?;
    values.resize(n, value);
    Ok(values)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn cell_error(index: usize, error: IonizationNetworkError) -> Error {
    if matches!(error, IonizationNetworkError::Cancelled) {
        Error::Cancelled
    } else {
        Error::Cell { index, error }
    }
}
