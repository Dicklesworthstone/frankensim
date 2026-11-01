//! Absolute sensor measurements in the native OpenEXR subset.
//!
//! FLOAT channels retain physical irradiance, illuminance, polarization, color
//! response and empirical sampling errors. There is no exposure, tone mapping,
//! denoising, peak normalization, or radiance-film conversion. Binary64 values
//! round to binary32; overflow and loss of a nonzero value to zero refuse.

use super::super::FiniteDetectorError;
use super::DetectorImage;
use fs_exec::Cx;
use fs_img::{Channel, ExrAttribute, ExrWriteLimits, ImgError, PixelType};

const CHANNELS: [&str; 20] = [
    "irradiance.I",
    "irradiance.Q",
    "irradiance.U",
    "irradiance.V",
    "illuminance.I",
    "illuminance.Q",
    "illuminance.U",
    "illuminance.V",
    "cie.X",
    "cie.Z",
    "stderr.irradiance.I",
    "stderr.irradiance.Q",
    "stderr.irradiance.U",
    "stderr.irradiance.V",
    "stderr.illuminance.I",
    "stderr.illuminance.Q",
    "stderr.illuminance.U",
    "stderr.illuminance.V",
    "stderr.cie.X",
    "stderr.cie.Z",
];

/// No partial artifact is returned on cancellation or a numerical refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum SensorImageExportError {
    Cancelled,
    Image(ImgError),
    Measurement(FiniteDetectorError),
    /// FLOAT storage cannot represent this finite physical channel value.
    UnrepresentableSample {
        pixel: usize,
        channel: &'static str,
    },
}
impl core::fmt::Display for SensorImageExportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sensor OpenEXR export: {self:?}")
    }
}
impl core::error::Error for SensorImageExportError {}
impl From<fs_exec::Cancelled> for SensorImageExportError {
    fn from(_: fs_exec::Cancelled) -> Self {
        Self::Cancelled
    }
}
impl From<ImgError> for SensorImageExportError {
    fn from(error: ImgError) -> Self {
        Self::Image(error)
    }
}
impl From<FiniteDetectorError> for SensorImageExportError {
    fn from(error: FiniteDetectorError) -> Self {
        Self::Measurement(error)
    }
}

impl DetectorImage {
    /// Encode physical pixel-area averages as an uncompressed FLOAT OpenEXR.
    ///
    /// `irradiance.*` and its errors are W/m²; `illuminance.*` are lux. CIE
    /// X/Z share the 683-scaled colorimetric convention in `xyz_flux`; CIE Y
    /// is `illuminance.I`. Stokes components keep the detector reference axes.
    /// Raster row zero remains the top row. Metadata retains cropped power,
    /// total detector power, quadrature, sample count and both physical pixel
    /// pitches. The native EXR display convention uses square pixels; the pitch
    /// metadata retains rectangular physical cells without changing their area.
    ///
    /// Writer limits are checked before allocating pixel planes. They bound
    /// encoded output and ordering references, not total process memory; the
    /// FLOAT projection planes require another 80 bytes per pixel. Cancellation
    /// is polled while projecting pixels and immediately before/after the
    /// native bounded encode, whose current API has no interior poll callback.
    pub fn to_openexr(
        &self,
        cx: &Cx<'_>,
        limits: ExrWriteLimits,
    ) -> Result<Vec<u8>, SensorImageExportError> {
        cx.checkpoint()?;
        let raster = self.raster();
        let count = raster.len();
        if self.pixels().len() != count {
            return Err(ImgError::Shape {
                expected: count,
                got: self.pixels().len(),
                context: "sensor pixel count",
            }
            .into());
        }
        let attributes = attributes(self);
        let layout = CHANNELS.map(|name| (name, PixelType::Float));
        let requirements = fs_img::exr_write_requirements_for_layout(
            raster.cols(),
            raster.rows(),
            &layout,
            &attributes,
        )?;
        if requirements.scratch_bytes > limits.max_scratch_bytes {
            return Err(ImgError::ResourceLimit {
                resource: "EXR ordering scratch",
                requested: requirements.scratch_bytes,
                limit: limits.max_scratch_bytes,
            }
            .into());
        }
        if requirements.output_bytes > limits.max_output_bytes {
            return Err(ImgError::ResourceLimit {
                resource: "EXR encoded output",
                requested: requirements.output_bytes,
                limit: limits.max_output_bytes,
            }
            .into());
        }
        let mut channels = Vec::new();
        channels
            .try_reserve_exact(CHANNELS.len())
            .map_err(|_| ImgError::AllocationRefused {
                resource: "sensor EXR channel descriptors",
                requested: (CHANNELS.len() * core::mem::size_of::<Channel>()) as u64,
            })?;
        for name in CHANNELS {
            cx.checkpoint()?;
            let mut data = Vec::new();
            data.try_reserve_exact(count)
                .map_err(|_| ImgError::AllocationRefused {
                    resource: "sensor EXR FLOAT plane",
                    requested: (count as u64) * 4,
                })?;
            channels.push(Channel {
                name: name.into(),
                ty: PixelType::Float,
                data,
            });
        }
        let area = raster.pixel_area_m2();
        for (index, pixel) in self.pixels().iter().enumerate() {
            if index % 256 == 0 {
                cx.checkpoint()?;
            }
            let average = pixel.area_average(area)?;
            let watts = average.radiant_w_m2;
            let lux = average.photopic_lux;
            let xyz = pixel.xyz_flux();
            let xyz_error = pixel.xyz_standard_error();
            let values = [
                watts.i,
                watts.q,
                watts.u,
                watts.v,
                lux.i,
                lux.q,
                lux.u,
                lux.v,
                xyz[0] / area,
                xyz[2] / area,
                average.radiant_standard_error_w_m2[0],
                average.radiant_standard_error_w_m2[1],
                average.radiant_standard_error_w_m2[2],
                average.radiant_standard_error_w_m2[3],
                average.photopic_standard_error_lux[0],
                average.photopic_standard_error_lux[1],
                average.photopic_standard_error_lux[2],
                average.photopic_standard_error_lux[3],
                xyz_error[0] / area,
                xyz_error[2] / area,
            ];
            for ((channel, value), name) in channels.iter_mut().zip(values).zip(CHANNELS) {
                let stored = value as f32;
                if !stored.is_finite() || (stored == 0.0 && value != 0.0) {
                    return Err(SensorImageExportError::UnrepresentableSample {
                        pixel: index,
                        channel: name,
                    });
                }
                channel.data.push(stored);
            }
        }
        cx.checkpoint()?;
        let bytes = fs_img::write_exr_with_attributes_budgeted(
            raster.cols(),
            raster.rows(),
            &channels,
            &attributes,
            limits,
        )?;
        cx.checkpoint()?;
        Ok(bytes)
    }
}

fn attributes(image: &DetectorImage) -> Vec<ExrAttribute> {
    let string = |name: &str, value: &str| ExrAttribute {
        name: name.into(),
        ty: "string".into(),
        value: value.as_bytes().to_vec(),
    };
    let double = |name: &str, value: f64| ExrAttribute {
        name: name.into(),
        ty: "double".into(),
        value: value.to_le_bytes().to_vec(),
    };
    vec![
        string(
            "frankensim.sensorScope",
            "geometric-optics finite-pixel area-average; raw Monte Carlo estimate; stderr is empirical, not a confidence bound",
        ),
        string(
            "frankensim.sensorUnits",
            "irradiance.*=W/m^2; illuminance.*=lux; cie.X/Z=683-scaled tristimulus per m^2; cie.Y=illuminance.I; stderr.*=matching channel units",
        ),
        string(
            "frankensim.sensorRows",
            "row0=top; columns=left-to-right; detector Stokes reference axes",
        ),
        string(
            "frankensim.sensorSpectrum",
            "80 midpoint bins, 380-780 nm, delta=5 nm; CIE 1931 2-degree observer; efficacy=683 lm/W",
        ),
        string(
            "frankensim.sensorSamples",
            &image.total_detector().sample_count().to_string(),
        ),
        double("frankensim.pixelArea_m2", image.raster().pixel_area_m2()),
        double(
            "frankensim.pixelWidth_m",
            image.raster().width_m() / f64::from(image.raster().cols()),
        ),
        double(
            "frankensim.pixelHeight_m",
            image.raster().height_m() / f64::from(image.raster().rows()),
        ),
        double(
            "frankensim.contactResolution_m",
            image.total_detector().max_contact_resolution_m(),
        ),
        double("frankensim.cropped_W", image.cropped().radiant_stokes_w().i),
        double(
            "frankensim.cropped_lm",
            image.cropped().photopic_stokes_lm().i,
        ),
        double(
            "frankensim.total_W",
            image.total_detector().radiant_stokes_w().i,
        ),
        double(
            "frankensim.total_lm",
            image.total_detector().photopic_stokes_lm().i,
        ),
    ]
}

#[cfg(test)]
mod tests;
