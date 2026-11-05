//! Calibrated point-source transmission through an ordinary native scene slab.
//!
//! The selected image crosses the two broad parallel glass faces exactly once.
//! The existing slab solver supplies the wavelength-specific Snell connection
//! and receiving area per source solid angle. Native s/p transmission preserves
//! the polarization induced at the first face through the second face; native
//! Beer absorption uses the actual internal path length. Source intensity is
//! evaluated along the solved launch direction, never the unrefracted chord.
//!
//! [`PhotometricSource::transmitted_slab_irradiance`] selects the image with
//! no internal returns. [`PhotometricSource::transmitted_slab_history_irradiance`]
//! sums the incoherent broad-face transmitted family, with explicit work caps
//! and an omitted-image bound. Neither includes other scene illumination,
//! side-face paths, diffuse/scattered light or coherent interference between
//! slab faces. The actual material's coherent surface coatings remain part of
//! each individual interface. No calibration or roundoff enclosure is inferred.

use super::bundle::{BundleIrradianceMeter, BundleIrradianceReading, PlanarIrradianceReceiver};
use super::{PhotometricSource, SPECTRAL_BINS, wavelength_nm};
use crate::gradient_index::bundle::ProjectionOptions;
use crate::photometry::{IlluminanceProbe, PhotometryError, TransportOptions};
use crate::tracer::{Scene, photometric::slab::SlabConnection};
use fs_exec::Cx;
use fs_geom::Vec3;

mod history;
pub use history::{SlabHistoryError, SlabHistoryOptions, TransmittedSlabHistorySpectrum};

/// Complete native 80-bin spectrum for the selected two-transmission image.
/// Opaque blockers can remove different bins when dispersion changes the path.
#[derive(Clone, Debug, PartialEq)]
pub struct TransmittedSlabSpectrum {
    bins_w_m2_nm: [f64; SPECTRAL_BINS],
    launch_directions: [Vec3; SPECTRAL_BINS],
    unoccluded: [bool; SPECTRAL_BINS],
    reading: BundleIrradianceReading,
}

impl TransmittedSlabSpectrum {
    /// Selected-image irradiance at the native bin midpoints, in W/(m² nm).
    pub const fn spectral_irradiance_w_m2_nm(&self) -> &[f64; SPECTRAL_BINS] {
        &self.bins_w_m2_nm
    }

    /// Actual source-to-slab launch direction for each native wavelength bin.
    pub const fn launch_directions(&self) -> &[Vec3; SPECTRAL_BINS] {
        &self.launch_directions
    }

    /// True only when all three physical segments reach the selected receiver.
    pub const fn unoccluded(&self) -> &[bool; SPECTRAL_BINS] {
        &self.unoccluded
    }

    /// Visible radiant irradiance and lux of this selected image only.
    pub const fn reading(&self) -> BundleIrradianceReading {
        self.reading
    }
}

impl PhotometricSource {
    /// Measure one directly transmitted image through a named ordinary scene
    /// primitive. The glass must be a closed static mesh instance with smooth,
    /// parallel broad faces admitted by the native slab connector. The source
    /// and receiver must lie in vacuum on opposite sides. Their straight chord
    /// must identify the selected broad face; the actual dispersive paths are
    /// subsequently solved and checked against every original scene primitive.
    ///
    /// `receiver.normal()` faces the arriving light. An opaque blocker gives a
    /// zero bin; rough, nested, nonparallel, edge-clipped, back-facing or
    /// unresolved connections refuse the complete reading. A receiver on an
    /// opaque surface uses the explicit endpoint exclusion in `options`.
    /// The distance budget bounds the full refracted path, not just its chord.
    ///
    /// Geometry is prepared once, then each wavelength uses fixed local
    /// storage with cancellation checks. Existing total path jobs keep their
    /// stricter point-source/glass admission: this selected-image measurement
    /// does not supply their missing reflected/specular source histories.
    pub fn transmitted_slab_irradiance(
        &self,
        cx: &Cx<'_>,
        scene: &Scene,
        slab_primitive: usize,
        receiver: IlluminanceProbe,
        options: TransportOptions,
    ) -> Result<TransmittedSlabSpectrum, PhotometryError> {
        let connection = SlabConnection::try_new(
            cx,
            scene,
            slab_primitive,
            self.position_m(),
            receiver,
            options,
        )?;
        let planar_receiver = PlanarIrradianceReceiver::try_new(
            receiver.position_m(),
            receiver.normal(),
            options.ray_epsilon_m,
        )
        .map_err(|_| PhotometryError::InvalidProbe)?;
        let mut meter =
            BundleIrradianceMeter::new(self, planar_receiver, ProjectionOptions::default());
        let mut result = TransmittedSlabSpectrum {
            bins_w_m2_nm: [0.0; SPECTRAL_BINS],
            launch_directions: [Vec3::new(0.0, 0.0, 0.0); SPECTRAL_BINS],
            unoccluded: [false; SPECTRAL_BINS],
            reading: BundleIrradianceReading {
                radiant_irradiance_w_m2: 0.0,
                illuminance_lux: 0.0,
            },
        };
        for bin in 0..SPECTRAL_BINS {
            cx.checkpoint()?;
            let wavelength = wavelength_nm(bin);
            let ray = connection.sample(cx, wavelength)?;
            let intensity = self.spectral_intensity_w_sr_nm(ray.launch_direction, wavelength)?;
            let value = ray.irradiance(intensity)?;
            meter
                .record_solved_irradiance(bin, value)
                .map_err(|_| PhotometryError::NonFiniteResult)?;
            result.bins_w_m2_nm[bin] = value;
            result.launch_directions[bin] = ray.launch_direction;
            result.unoccluded[bin] = ray.unoccluded;
        }
        result.reading = meter
            .reading()
            .map_err(|_| PhotometryError::NonFiniteResult)?;
        cx.checkpoint()?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
