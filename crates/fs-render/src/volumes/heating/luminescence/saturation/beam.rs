//! Depleting spectral pumps on the original cell and thermal timeline.
//!
//! Mutually incoherent, quadrature-weighted pumps share one finite inventory.
//! Each beam is transported in its own native DDA order.
//! Its molecular optical depth is cross section times the complete cell's
//! ground inventory times ray length / cell volume. The coupled thermal owner
//! solves that inventory and the transmitted photons together. This is an
//! explicitly cell-uniform finite-volume closure: strong longitudinal or
//! transverse saturation requires spatial refinement. No beam area is invented.
//! Molecular cross sections and kinetic rates refresh at the native thermal
//! interval; passive host coefficients retain the original optical snapshot.
//! Inner refinement controls population/opacity coupling, not temperature
//! feedback error.

use super::{
    CheckedSum, ContentHash, DepositionError, DepositionSimulation, DomainHasher,
    OpticalHeatingError, SaturableLuminescentHeating, finite, invalid, nonnegative,
    preserve_checkpoint, resolve_model,
};
use crate::volumes::heating::feedback::IncidentBeam;
use crate::volumes::transport::{AbsorbingGrid, VolumeTransportError, positive_product};
use fs_exec::Cx;
use fs_thermal::enthalpy::CellId;
use fs_thermal::luminescence::transfer::pump_beams::{
    SaturablePumpBeams, SaturablePumpPath, SaturablePumpSegment,
};
use fs_thermal::luminescence::transfer::pump_column::{SaturablePumpCell, SaturablePumpColumn};
use fs_thermal::luminescence::transfer::saturation::SaturableSingletTripletStep;
use std::sync::Arc;

pub use fs_thermal::luminescence::transfer::pump_beams::PumpBeamsLimits;
pub use fs_thermal::luminescence::transfer::pump_column::PumpColumnLimits;

mod accounts;
mod binding;
mod evolution;
pub(super) use accounts::BeamAccounts;
pub use accounts::{DepletingPumpBeamEnergy, DepletingPumpEnergy, DepletingPumpReport};
pub use binding::DepletingPumpBeam;

use super::cross_sections::{CrossSections, EmitterAbsorptionCrossSections};

/// The two incident pump nodes use the same molecular source carrier as
/// reabsorbed luminescence. Bulk mass absorption remains a passive host channel.
pub type BeamEmitterCrossSections = EmitterAbsorptionCrossSections<2>;

#[derive(Clone, Debug)]
struct BeamCell {
    /// Index in the unique shared-inventory vector, not this beam's path index.
    cell_index: usize,
    length_m: f64,
    volume_m3: f64,
    passive_optical_depth: [f64; 2],
    cross_sections: BeamEmitterCrossSections,
}

#[derive(Debug)]
struct BeamPath {
    source_identity: ContentHash,
    wavelengths_nm: [f64; 2],
    photon_energies_j: [f64; 2],
    incident_photon_rates_s: [f64; 2],
    pump_window_s: Option<(f64, f64)>,
    cells: Vec<BeamCell>,
}

#[derive(Debug)]
struct BeamModel {
    /// Original owner emitter indices, unique and in native cell order.
    cells: Vec<usize>,
    paths: Vec<BeamPath>,
    limits: PumpBeamsLimits,
}

#[derive(Clone, Debug)]
pub(super) struct DepletingBeam {
    model: Arc<BeamModel>,
    prepared: Option<Arc<PreparedPumps>>,
    accounts: BeamAccounts,
}

#[derive(Debug)]
enum PreparedSolver {
    /// Preserve the original one-column algorithm and work for one beam.
    Column(SaturablePumpColumn<2>),
    Beams(SaturablePumpBeams<2>),
}

/// Only immutable constitutive data is shared across accepted snapshots.
/// Populations, pulse activity and photon histories are never cached here.
#[derive(Debug)]
pub(super) struct PreparedPumps {
    temperatures: Vec<u64>,
    solver: PreparedSolver,
    /// Same original paths, retained for joint emitted-photon feedback.
    paths: Vec<SaturablePumpPath<2>>,
}

pub(super) struct BeamCandidate {
    pub steps: Vec<Option<SaturableSingletTripletStep>>,
    pub passive_heat_j: Vec<(CellId, f64)>,
    pub absolute_roundoff_j: f64,
    pub accounts: BeamAccounts,
}
