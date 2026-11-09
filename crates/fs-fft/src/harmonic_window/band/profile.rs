//! Positive factors for the actual normalized detector profile.
//!
//! Uniform and triangular profiles share source admission, DD phase clocks,
//! fixed storage and all contractions. Only the orthogonal temporal basis
//! changes. Existing uniform-only constructors retain their original domain.

use super::basis::{self, Units};
use super::*;

/// Owned positive factor for a uniform or triangular normalized observation.
/// The storage and contractions are shared with the original uniform factor.
pub type PreparedHarmonicBand = PreparedUniformHarmonicBand;
/// Fixed reusable storage for either admitted detector profile.
pub type HarmonicBandWorkspace = UniformHarmonicBandWorkspace;
/// Borrowed positive factor; mode zero retains the weighted actual-real mean.
pub type HarmonicBandRef<'a> = UniformHarmonicBandRef<'a>;

impl HarmonicWindow {
    /// Analytic-factor admission on original Hz for this detector profile.
    /// Admits 2..64 tones and a half-band phase at most 0.5; at most 17 modes.
    /// A common carrier cancels only from second moments. Independent source
    /// banks remain separate, and proper peak moments require scale 0.5.
    pub fn band_mode_count_hz(self, frequencies: &[f64]) -> Result<Option<usize>, WindowError> {
        Ok(basis::admit_profile(self, frequencies, Units::Hertz, false)?.map(|s| s.modes))
    }

    /// Analytic-factor admission on the original angular-frequency clock.
    pub fn band_mode_count_rad_s(self, frequencies: &[f64]) -> Result<Option<usize>, WindowError> {
        Ok(basis::admit_profile(self, frequencies, Units::Radians, false)?.map(|s| s.modes))
    }

    /// Actual-real admission: 2..32 tones and max(abs(2*pi*f))*T/2 <= 0.5.
    /// Absolute carrier and sum-frequency contributions remain present.
    pub fn real_band_mode_count_hz(
        self,
        frequencies: &[f64],
    ) -> Result<Option<usize>, WindowError> {
        Ok(basis::admit_profile(self, frequencies, Units::Hertz, true)?.map(|s| s.modes))
    }

    /// Actual-real admission with max(abs(omega))*T/2 <= 0.5.
    pub fn real_band_mode_count_rad_s(
        self,
        frequencies: &[f64],
    ) -> Result<Option<usize>, WindowError> {
        Ok(basis::admit_profile(self, frequencies, Units::Radians, true)?.map(|s| s.modes))
    }

    /// Cache an analytic positive factor with this window's actual profile.
    pub fn prepare_band_hz(
        self,
        frequencies: &[f64],
    ) -> Result<Option<PreparedHarmonicBand>, WindowError> {
        basis::prepare_profile_owned(self, frequencies, Units::Hertz, false)
    }

    /// Cache an analytic factor without converting original rad/s through Hz.
    pub fn prepare_band_rad_s(
        self,
        frequencies: &[f64],
    ) -> Result<Option<PreparedHarmonicBand>, WindowError> {
        basis::prepare_profile_owned(self, frequencies, Units::Radians, false)
    }

    /// Cache actual-real products and weighted mean on the original Hz clock.
    pub fn prepare_real_band_hz(
        self,
        frequencies: &[f64],
    ) -> Result<Option<PreparedHarmonicBand>, WindowError> {
        basis::prepare_profile_owned(self, frequencies, Units::Hertz, true)
    }

    /// Cache actual-real products and weighted mean on the original rad/s clock.
    pub fn prepare_real_band_rad_s(
        self,
        frequencies: &[f64],
    ) -> Result<Option<PreparedHarmonicBand>, WindowError> {
        basis::prepare_profile_owned(self, frequencies, Units::Radians, true)
    }

    /// Prepare this profile's analytic Hz factor without heap allocation.
    pub fn prepare_band_hz_in<'a>(
        self,
        frequencies: &[f64],
        workspace: &'a mut HarmonicBandWorkspace,
    ) -> Result<Option<HarmonicBandRef<'a>>, WindowError> {
        workspace::prepare_profile_in(self, frequencies, Units::Hertz, false, workspace)
    }

    /// Prepare this profile's analytic rad/s factor without heap allocation.
    pub fn prepare_band_rad_s_in<'a>(
        self,
        frequencies: &[f64],
        workspace: &'a mut HarmonicBandWorkspace,
    ) -> Result<Option<HarmonicBandRef<'a>>, WindowError> {
        workspace::prepare_profile_in(self, frequencies, Units::Radians, false, workspace)
    }

    /// Prepare actual-real Hz products with this profile in fixed scratch.
    pub fn prepare_real_band_hz_in<'a>(
        self,
        frequencies: &[f64],
        workspace: &'a mut HarmonicBandWorkspace,
    ) -> Result<Option<HarmonicBandRef<'a>>, WindowError> {
        workspace::prepare_profile_in(self, frequencies, Units::Hertz, true, workspace)
    }

    /// Prepare actual-real rad/s products with this profile in fixed scratch.
    pub fn prepare_real_band_rad_s_in<'a>(
        self,
        frequencies: &[f64],
        workspace: &'a mut HarmonicBandWorkspace,
    ) -> Result<Option<HarmonicBandRef<'a>>, WindowError> {
        workspace::prepare_profile_in(self, frequencies, Units::Radians, true, workspace)
    }
}
