//! Frequency-resolved observations of the original positive-frequency fields.

use super::*;

/// Borrowed complex analyzer weights at explicitly authored angular frequencies.
///
/// Each entry is `(omega_rad_s, projection)`. Frequencies must be finite,
/// positive, strictly increasing and unique. A query requires an EXACT entry
/// for every retained frequency in its original mean and factor banks; extra
/// entries are harmless. There is no interpolation or extrapolation of complex
/// phase, polarization angle, gain or material response.
/// Each response acts on all spatial modes in the original Cartesian basis;
/// angle-dependent optical scattering remains with the physical field owners.
///
/// The observed signals are `Re(sum_j a(omega_j)† E_j exp(i omega_j t))`
/// and `Re(sum_j b(omega_j)† H_j exp(i omega_j t))`. Absolute spectral phase is
/// retained, so dispersive retardance and common group delay affect finite
/// shutters and beats. For example, weights proportional to `exp(i omega tau)`
/// observe the original signal delayed by `tau` under the native +i omega t
/// convention. The dagger conjugates the supplied weights.
///
/// The caller supplies the actual calibration/material response at those
/// frequencies. Finite line samples do not establish passivity, causality,
/// out-of-band behavior, detector absorption or source back-action. Weights
/// retain the constant projection's dimensions and normalization convention.
/// Admission is cancellable and allocates nothing; the caller owns the table.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumExposureApertureJonesSpectrum<'a> {
    samples: &'a [(f64, AngularSpectrumExposureApertureJonesProjection)],
}

impl<'a> AngularSpectrumExposureApertureJonesSpectrum<'a> {
    pub fn try_new(
        samples: &'a [(f64, AngularSpectrumExposureApertureJonesProjection)],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if samples.is_empty() {
            return Err(bad("spectral analyzer needs at least one frequency"));
        }
        let mut previous = 0.0;
        for (index, &(omega, _)) in samples.iter().enumerate() {
            if index % 32 == 0 {
                poll(&mut cancelled)?;
            }
            if !omega.is_finite() || omega <= previous {
                return Err(bad(
                    "spectral analyzer frequencies must be finite, positive and strictly increasing",
                ));
            }
            previous = omega;
        }
        poll(&mut cancelled)?;
        Ok(Self { samples })
    }

    /// The original frequency/response pairs, without resampling or allocation.
    pub fn samples(self) -> &'a [(f64, AngularSpectrumExposureApertureJonesProjection)] {
        self.samples
    }

    /// Exact-frequency response. A missing frequency refuses; nearby samples
    /// do not supply an implicit phase interpolation or material continuation.
    pub fn projection_at(
        self,
        angular_frequency_rad_s: f64,
    ) -> Result<AngularSpectrumExposureApertureJonesProjection, AngularSpectrumError> {
        self.samples
            .binary_search_by(|(omega, _)| omega.total_cmp(&angular_frequency_rad_s))
            .map(|index| self.samples[index].1)
            .map_err(|_| bad("spectral analyzer does not cover a source frequency"))
    }
}

impl ProjectionLaw for AngularSpectrumExposureApertureJonesSpectrum<'_> {
    fn admit_frequencies(
        &self,
        angular_frequencies_rad_s: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AngularSpectrumError> {
        for (index, &omega) in angular_frequencies_rad_s.iter().enumerate() {
            if index % 32 == 0 {
                poll(cancelled)?;
            }
            self.projection_at(omega)?;
        }
        poll(cancelled)
    }

    fn at(
        &self,
        angular_frequency_rad_s: f64,
    ) -> Result<(AngularSpectrumExposureApertureJonesProjection, bool), AngularSpectrumError> {
        let projection = self.projection_at(angular_frequency_rad_s)?;
        Ok((projection, real_weights(projection)))
    }
}

impl<'a> Projected<AngularSpectrumExposureApertureJonesSpectrum<'a>> {
    pub(super) fn spectral(projection: AngularSpectrumExposureApertureJonesSpectrum<'a>) -> Self {
        Self {
            projection,
            sums: [Dd::ZERO; 3],
        }
    }
}

impl AngularSpectrumExposure<'_, '_, '_> {
    /// Apply an authored complex frequency response to each original grouped
    /// DD field before products. Complete frequency coverage is checked before
    /// spatial caches or Maxwell reconstruction. Carrier sums and beats remain.
    #[allow(clippy::too_many_arguments)]
    pub fn measure_aperture_spectral_projection(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        projection: AngularSpectrumExposureApertureJonesSpectrum<'_>,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureProjectionReading, AngularSpectrumError> {
        measure(
            self.anchor,
            || core::iter::once((self, false)),
            ExposureSide::OneSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            Projected::spectral(projection),
            &mut cancelled,
        )
    }

    /// Frequency-dependent field analyzers at either original illumination
    /// face, retaining its complete incoming/outgoing phase relation.
    #[allow(clippy::too_many_arguments)]
    pub fn measure_illumination_aperture_spectral_projection(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        projection: AngularSpectrumExposureApertureJonesSpectrum<'_>,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureProjectionReading, AngularSpectrumError> {
        measure(
            self.anchor,
            || core::iter::once((self, false)),
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            Projected::spectral(projection),
            &mut cancelled,
        )
    }
}

impl AngularSpectrumExposureEnsemble<'_, '_, '_, '_> {
    /// The same authored spectral analyzer acts on the mean and every real or
    /// proper-complex factor. Every bank's retained frequencies must be covered;
    /// independent factors retain their original statistics and add as moments.
    #[allow(clippy::too_many_arguments)]
    pub fn measure_aperture_spectral_projection(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        projection: AngularSpectrumExposureApertureJonesSpectrum<'_>,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureProjectionReading, AngularSpectrumError> {
        measure(
            self.anchor,
            || self.banks(),
            ExposureSide::OneSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            Projected::spectral(projection),
            &mut cancelled,
        )
    }

    /// Expected frequency-dependent analyzer moments at either original face.
    #[allow(clippy::too_many_arguments)]
    pub fn measure_illumination_aperture_spectral_projection(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        projection: AngularSpectrumExposureApertureJonesSpectrum<'_>,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureProjectionReading, AngularSpectrumError> {
        measure(
            self.anchor,
            || self.banks(),
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            Projected::spectral(projection),
            &mut cancelled,
        )
    }
}
