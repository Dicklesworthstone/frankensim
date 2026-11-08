//! Visible observation of wavelength-resolved direct-beam opacity.

use super::{DvrError, MAX_OPACITY_WAVELENGTHS, SpectralOpacityImage, reserve_image};
use crate::spectral::{
    observer::{XyzAccumulator, XyzObserverError, XyzObserverLane},
    y_integral,
};

impl SpectralOpacityImage {
    /// Approximate visible XYZ of a supplied spectrum transmitted by the image.
    ///
    /// Background samples are spectral radiance per nanometre at the image's
    /// original wavelengths; `quadrature_weights_nm` supplies the corresponding
    /// nonnegative integration widths in nanometres. No illuminant, bin width,
    /// spectral interpolation, exposure, or display transform is inferred.
    /// The same spectrum illuminates every pixel and is constant within each
    /// pixel footprint. Thus a cold, colored absorber is observed through its
    /// transmitted light instead of being assigned an emission color.
    ///
    /// This uses the renderer's analytic CIE fits, limited to 380--780 nm, and
    /// its relative `1 / integral(y_bar)` normalization. It is an approximate
    /// color observation, not calibrated luminance or a broadband power result.
    /// The original spectral samples remain available, including infrared.
    pub fn approximate_transmitted_xyz(
        &self,
        quadrature_weights_nm: &[f64],
        background_w_m2_sr_nm: &[f64],
    ) -> Result<Vec<[f64; 3]>, DvrError> {
        if background_w_m2_sr_nm.len() != self.lanes {
            return Err(DvrError::InvalidBackgroundSize);
        }
        self.transmitted_xyz(quadrature_weights_nm, |_| background_w_m2_sr_nm)
    }

    /// Approximate visible XYZ against one supplied background spectrum per pixel.
    ///
    /// Background layout is pixel-major, wavelength-fastest, matching
    /// [`Self::samples`]. It must contain exactly one radiance-per-nanometre
    /// sample per stored opacity sample. Each supplied spectrum must be constant
    /// over its own pixel footprint: unresolved background/occluder correlation
    /// cannot be recovered from averaged transmission. Units, observer limits,
    /// and normalization match [`Self::approximate_transmitted_xyz`].
    ///
    /// Projection writes XYZ directly, without another full spectral image.
    /// Retained optical depth preserves bright-background contributions even
    /// when unit transmittance has underflowed. Any invalid input or arithmetic
    /// refusal returns no partial image.
    pub fn approximate_transmitted_xyz_over_image(
        &self,
        quadrature_weights_nm: &[f64],
        backgrounds_w_m2_sr_nm: &[f64],
    ) -> Result<Vec<[f64; 3]>, DvrError> {
        if backgrounds_w_m2_sr_nm.len() != self.samples.len() {
            return Err(DvrError::InvalidBackgroundSize);
        }
        self.transmitted_xyz(quadrature_weights_nm, |pixel| {
            let first = pixel * self.lanes;
            &backgrounds_w_m2_sr_nm[first..first + self.lanes]
        })
    }

    fn transmitted_xyz<'a>(
        &self,
        quadrature_weights_nm: &[f64],
        background: impl Fn(usize) -> &'a [f64],
    ) -> Result<Vec<[f64; 3]>, DvrError> {
        if quadrature_weights_nm.len() != self.lanes {
            return Err(DvrError::InvalidSpectralWeights);
        }
        let normalization = 1.0 / y_integral();
        let mut observer = [XyzObserverLane::ZERO; MAX_OPACITY_WAVELENGTHS];
        for lane in 0..self.lanes {
            observer[lane] = XyzObserverLane::new(
                self.wavelengths_nm[lane],
                quadrature_weights_nm[lane],
                normalization,
            )
            .map_err(|error| observer_error(error, lane))?;
        }
        let mut image = reserve_image(self.samples.len() / self.lanes)?;
        for (index, pixel) in self.samples.chunks_exact(self.lanes).enumerate() {
            let background = background(index);
            let mut xyz = XyzAccumulator::default();
            for lane in 0..self.lanes {
                xyz.add(
                    &observer[lane],
                    0.0,
                    pixel[lane].transmittance,
                    pixel[lane].optical_depth,
                    background[lane],
                )
                .map_err(|error| observer_error(error, lane))?;
            }
            image.push(xyz.finish().map_err(|error| observer_error(error, 0))?);
        }
        Ok(image)
    }
}

fn observer_error(error: XyzObserverError, lane: usize) -> DvrError {
    match error {
        XyzObserverError::Wavelength => DvrError::InvalidObserverWavelength { lane },
        XyzObserverError::Weight => DvrError::InvalidSpectralWeights,
        XyzObserverError::Radiance => DvrError::InvalidBackgroundRadiance,
        XyzObserverError::Unrepresentable => DvrError::UnrepresentableRadiance,
    }
}

#[cfg(test)]
mod tests;
