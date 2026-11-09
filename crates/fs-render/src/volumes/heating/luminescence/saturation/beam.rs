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
use fs_thermal::luminescence::transfer::pump_column::{
    SaturablePumpCell, SaturablePumpColumn,
};
use fs_thermal::luminescence::transfer::saturation::SaturableSingletTripletStep;
use std::sync::Arc;
use fs_thermal::luminescence::transfer::pump_beams::{
    SaturablePumpBeams, SaturablePumpPath, SaturablePumpSegment,
};

pub use fs_thermal::luminescence::transfer::pump_column::PumpColumnLimits;
pub use fs_thermal::luminescence::transfer::pump_beams::PumpBeamsLimits;

mod accounts;
mod binding;
mod evolution;
pub use accounts::{DepletingPumpBeamEnergy, DepletingPumpEnergy, DepletingPumpReport};
pub use binding::DepletingPumpBeam;
pub(super) use accounts::BeamAccounts;

#[derive(Clone, Debug)]
enum CrossSections {
    Authored {
        identity: ContentHash,
        values: [[f64; 2]; 2],
    },
    #[cfg(feature = "material-state")]
    Material {
        card: Arc<fs_matdb::MaterialCard>,
        source: fs_material::phosphorescence::pump::GroundStatePumpStatePoint,
    },
}

/// Molecular cross sections for one original emitter cell. Rows are the two
/// incident wavelength lanes; columns are transitions [ground->singlet,
/// ground->triplet]. Both transitions compete for the SAME incident photons.
/// Bulk mass absorption is a separate passive host channel.
#[derive(Clone, Debug)]
pub struct BeamEmitterCrossSections {
    grid_cell: usize,
    source: CrossSections,
}

impl BeamEmitterCrossSections {
    /// Bind explicit, nonnegative square metres per ground-state emitter.
    /// The source identity attributes authored data; it does not certify it.
    pub fn authored(
        grid_cell: usize,
        source_identity: ContentHash,
        cross_sections_m2: [[f64; 2]; 2],
    ) -> Result<Self, OpticalHeatingError> {
        if source_identity == ContentHash([0; 32])
            || !cross_sections_m2.into_iter().flatten().all(nonnegative)
        {
            return Err(invalid("invalid molecular pump cross sections"));
        }
        Ok(Self {
            grid_cell,
            source: CrossSections::Authored {
                identity: source_identity,
                values: cross_sections_m2,
            },
        })
    }

    /// Bind the original selected molecular cross-section claims. Both vacuum
    /// wavelength and absolute temperature must be explicit owned coordinates;
    /// every other physical condition remains fixed by the native query.
    #[cfg(feature = "material-state")]
    pub fn from_material(
        grid_cell: usize,
        card: Arc<fs_matdb::MaterialCard>,
        source: fs_material::phosphorescence::pump::GroundStatePumpStatePoint,
    ) -> Result<Self, OpticalHeatingError> {
        if source.temperature_axis().is_none() {
            return Err(invalid(
                "thermal pump source requires an owned absolute-temperature axis",
            ));
        }
        // Re-query at the exact initial point to reject a different card now.
        source
            .at_query_point(&card, &source.resolved().query())
            .map_err(|_| invalid("pump cross sections differ from original material card"))?;
        Ok(Self {
            grid_cell,
            source: CrossSections::Material { card, source },
        })
    }

    pub const fn grid_cell(&self) -> usize {
        self.grid_cell
    }

    fn identity(&self) -> ContentHash {
        match &self.source {
            CrossSections::Authored { identity, .. } => *identity,
            #[cfg(feature = "material-state")]
            CrossSections::Material { source, .. } => source.identity(),
        }
    }

    fn evaluate(
        &self,
        wavelengths_nm: [f64; 2],
        temperature_k: f64,
    ) -> Result<[[f64; 2]; 2], OpticalHeatingError> {
        match &self.source {
            CrossSections::Authored { values, .. } => {
                let _ = (wavelengths_nm, temperature_k);
                Ok(*values)
            }
            #[cfg(feature = "material-state")]
            CrossSections::Material { card, source } => {
                let mut values = [[0.0; 2]; 2];
                for lane in 0..2 {
                    values[lane] = source.at_wavelength_temperature(
                        card,
                        finite(wavelengths_nm[lane] / 1e9)?,
                        Some(temperature_k),
                    ).map_err(|_| invalid("molecular pump wavelength/temperature is outside original material support"))?
                        .cross_sections_m2();
                }
                Ok(values)
            }
        }
    }
}

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
}

pub(super) struct BeamCandidate {
    pub steps: Vec<Option<SaturableSingletTripletStep>>,
    pub passive_heat_j: Vec<(CellId, f64)>,
    pub absolute_roundoff_j: f64,
    pub accounts: BeamAccounts,
}
