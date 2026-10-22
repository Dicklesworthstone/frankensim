//! Material walls and gas fields share one regional acoustic binding.
//!
//! The material owner supplies the original thin-membrane reduction and
//! generalized-Maxwell spectrum. fs-duct supplies local geometry, gas loss and
//! propagation. Feedback and outlet observation cannot select different wall
//! parameters from this binding.

use super::{RegionalDuct, DrivingPointError, LineRealization};
use super::outlet::{OutletError, OutletPressure};
use crate::acoustic_realize::AcousticRealizeError;
use crate::pcm_wav::baffled::CircularOutletReceiver;
pub use crate::bernoulli_aperture::wall::material::LocalMaterialWall;
use fs_duct::{DuctError, DuctResponse, LossModel, Termination, WallLaw, WallModel};
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
    pub fn characteristic_line(&self, termination: Termination, sample_rate_hz: u32,
        n: usize, n_fft_override: Option<usize>) -> Result<DelayedFilter, DrivingPointError> {
        self.duct.characteristic_line_with_wall(termination, sample_rate_hz, n,
            WallModel::Law(&self.wall), n_fft_override)
    }

    /// Same feedback realization, retaining its original discretization receipt.
    pub fn characteristic_line_realized(&self, termination: Termination, sample_rate_hz: u32,
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
