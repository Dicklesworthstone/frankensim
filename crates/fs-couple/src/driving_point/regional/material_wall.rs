//! Material walls and gas fields share one regional acoustic binding.
//!
//! The material owner supplies the original thin-membrane reduction and
//! generalized-Maxwell spectrum. fs-duct supplies local geometry, gas loss and
//! propagation. Feedback and outlet observation cannot select different wall
//! parameters from this binding.

use super::{RegionalDuct, DrivingPointError, LineRealization, Region, invalid};
use super::outlet::{OutletError, OutletPressure};
use crate::acoustic_realize::AcousticRealizeError;
use crate::pcm_wav::baffled::CircularOutletReceiver;
pub use crate::bernoulli_aperture::wall::material::LocalMaterialWall;
use fs_duct::{Duct, DuctError, DuctResponse, LossModel, MouthLoad, Segment, Termination, WallLaw, WallModel};
use fs_math::c64::C64;
use fs_vfit::discretize::DelayedFilter;

impl WallLaw for LocalMaterialWall {
    fn specific_impedance(&self, radius_m: f64, omega: f64) -> Result<C64, DuctError> {
        LocalMaterialWall::specific_impedance(self, radius_m, omega).map_err(|error| match error {
            AcousticRealizeError::InvalidDescription { what } => DuctError::BadParameter { what },
            _ => DuctError::BadParameter { what: "source-bound material wall impedance is not representable" },
        })
    }
}

/// Complete source-owned wall assignment, independent of the gas partition.
/// Each segment keeps its own material state, thickness, restraint and memory.
/// No effective modulus or material-name inference is introduced at a join.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialWallLayout {
    /// One original local membrane on the complete admitted duct.
    Uniform(LocalMaterialWall),
    /// Exactly one source per authored axial segment, in geometry order.
    /// This bounded image excludes tone holes and continuous gas-field bindings.
    Segments(Vec<LocalMaterialWall>),
}

impl MaterialWallLayout {
    /// Uniform source, when that assignment was explicitly selected.
    #[must_use]
    pub const fn uniform(&self) -> Option<&LocalMaterialWall> {
        match self { Self::Uniform(wall) => Some(wall), Self::Segments(_) => None }
    }

    /// Original complete segment partition, without gas-run coalescing.
    #[must_use]
    pub fn segments(&self) -> Option<&[LocalMaterialWall]> {
        match self { Self::Segments(walls) => Some(walls), Self::Uniform(_) => None }
    }
}

/// One immutable gas/geometry/material binding for all acoustic consumers.
///
/// Geometry is admitted by the original regional constructor. Every actual
/// local radius/frequency is checked by the source wall during evaluation.
/// The same supplied wall material and thickness cover the axial duct and any
/// admitted chimneys. Different wall regions and axial shell coupling are not
/// inferred. Existing nonuniform tone-hole restrictions remain unchanged.
#[derive(Clone, Debug)]
pub struct MaterialWallDuct {
    duct: RegionalDuct,
    wall: LocalMaterialWall,
}

impl RegionalDuct {
    /// Attach an already source-resolved local membrane law to every consumer.
    /// This cold binding does not reset or mutate any running physical state.
    #[must_use]
    pub fn with_material_wall(self, wall: LocalMaterialWall) -> MaterialWallDuct {
        MaterialWallDuct { duct: self, wall }
    }

    /// Bind a complete axial material partition to the original regional TMM.
    /// Adjacent equal gas AND wall states coalesce. Every other interface keeps
    /// the existing pressure/volume-flow continuity, with no averaged material.
    /// Feedback, arbitrary terminal loads and outlet observations then use the
    /// retained walls when passed `WallModel::Rigid`/no external wall override.
    /// The rigid argument means no additional wall law, not a rigid substitute.
    ///
    /// This is cold, fixed-geometry, locally reacting acoustics. Wall joins add
    /// no axial shell/joint dynamics. Continuous gas fields and tone holes need
    /// a spatial wall/mutual-mouth owner and explicitly refuse on this path.
    pub fn with_material_wall_regions(mut self, walls: &[LocalMaterialWall])
        -> Result<Self, DuctError> {
        if self.continuous.is_some() || self.has_material_wall_regions() {
            return Err(invalid("segment wall binding cannot replace a continuous gas field or an existing wall partition"));
        }
        if walls.len() != self.gases.len() {
            return Err(invalid("material wall partition needs exactly one source per axial segment"));
        }
        let mut regions: Vec<Region> = Vec::new();
        for ((segment, gas), wall) in self.regions.iter().flat_map(|r| &r.duct.segments)
            .zip(&self.gases).zip(walls) {
            if matches!(segment, Segment::ToneHole { .. }) {
                return Err(invalid("segment material walls do not supply a cross-region tone-hole law"));
            }
            for radius in [segment.inlet_radius(), segment.outlet_radius()] {
                wall.base_patch(radius, segment.axial_length()).map_err(|_| {
                    invalid("segment geometry lies outside its material wall support")
                })?;
            }
            if regions.last().is_none_or(|r| r.gas != *gas || r.wall.as_ref() != Some(wall)) {
                regions.push(Region { duct: Duct { segments: Vec::new() },
                    gas: *gas, wall: Some(wall.clone()) });
            }
            regions.last_mut().expect("admitted material run").duct.segments.push(segment.clone());
        }
        self.regions = regions;
        Ok(self)
    }

    pub(super) fn has_material_wall_regions(&self) -> bool {
        self.regions[0].wall.is_some()
    }

    pub(super) fn admit_wall_override(&self, wall: WallModel<'_>) -> Result<(), DuctError> {
        if self.has_material_wall_regions() && !wall.is_rigid() {
            return Err(invalid("a retained material wall partition cannot also select a wall override"));
        }
        Ok(())
    }

    pub(super) fn region_wall<'a>(region: &'a Region, fallback: WallModel<'a>) -> WallModel<'a> {
        match &region.wall { Some(wall) => WallModel::Law(wall), None => fallback }
    }
}

impl MaterialWallDuct {
    /// Original media, geometry, exterior and field-resolution choices.
    #[must_use] pub const fn duct(&self) -> &RegionalDuct { &self.duct }
    /// Same source identity, geometry reduction and spectrum used by every port.
    #[must_use] pub const fn wall(&self) -> &LocalMaterialWall { &self.wall }

    /// Complete input impedance and physical outlet transfer.
    pub fn response(&self, omega: f64, loss: LossModel, termination: Termination)
        -> Result<DuctResponse, DuctError> {
        self.duct.response_with_wall(omega, loss, termination, WallModel::Law(&self.wall))
    }

    /// Direct inlet reflection, retaining the continuous/stable propagation path.
    pub fn reflection(&self, omega: f64, loss: LossModel, termination: Termination)
        -> Result<C64, DuctError> {
        self.duct.reflection_with_wall(omega, loss, termination, WallModel::Law(&self.wall))
    }

    /// The existing finite-history characteristic feedback realization.
    /// Every sampled frequency must lie in the original material's support.
    /// Maxwell midpoint-state limits belong to the graph integrator; this
    /// frequency-domain image uses the original FIR sampling/history controls.
    /// A tabulated mouth also requires wall support at the existing 1e-3 rad/s
    /// quasi-static bore probe. Positive FFT-band coverage alone does not
    /// authorize unavailable static wall data; that case explicitly refuses.
    pub fn characteristic_line<'a>(&self, termination: impl Into<MouthLoad<'a>>, sample_rate_hz: u32,
        n: usize, n_fft_override: Option<usize>) -> Result<DelayedFilter, DrivingPointError> {
        self.duct.characteristic_line_with_wall(termination, sample_rate_hz, n,
            WallModel::Law(&self.wall), n_fft_override)
    }

    /// Same feedback realization, retaining its original discretization receipt.
    pub fn characteristic_line_realized<'a>(&self, termination: impl Into<MouthLoad<'a>>, sample_rate_hz: u32,
        n: usize, n_fft_override: Option<usize>)
        -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
        self.duct.characteristic_line_with_wall_realized(termination, sample_rate_hz, n,
            WallModel::Law(&self.wall), n_fft_override)
    }

    /// Actual baffled outlet pressure using the identical wall and exterior.
    /// This observes the outlet; it does not invent external wall radiation.
    pub fn baffled_outlet(&self, sample_rate_hz: u32, fft_samples: usize,
        receiver: CircularOutletReceiver) -> Result<OutletPressure, OutletError> {
        self.duct.baffled_outlet_with_wall(sample_rate_hz, fft_samples, receiver,
            WallModel::Law(&self.wall))
    }
}

#[cfg(test)]
mod tests;
