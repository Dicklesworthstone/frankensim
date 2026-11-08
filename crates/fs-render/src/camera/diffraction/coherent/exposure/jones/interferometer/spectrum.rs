//! Frequency-indexed passive arm optics on the original shared source spectrum.
//!
//! Preparation resolves each exact authored frequency once. Observation indexes
//! the retained sequence directly, preserving relative phase, filtering and
//! retardance before the original spectral/polarized interference contraction.

use super::super::{CoherentDiffractionExposurePixelLimits, CoherentDiffractionPixel};
use super::{
    CoherentDiffractionError, CoherentDiffractionInterferometer,
    CoherentInterferometerExposureLimits, CoherentJonesArm, InterferometerArms,
    PolarizedCoherentInterferometerExposureSample, PolarizedCoherentInterferometerPixelExposure,
    PreparedCoherentJonesDiffractionExposure, admit_count, limit, poll,
};

/// Explicit passive arm response at one original vacuum frequency.
///
/// Every distinct source frequency requires exactly one response, shared by all
/// source entries at that frequency. The arm sequences may have frequency-specific
/// losses, retardance, axes and common phases. Each sequence remains spatially
/// uniform; there is no interpolation, extrapolation or inferred material model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentJonesInterferometerResponse {
    pub frequency_hz: f64,
    pub arms: [CoherentJonesArm; 2],
}

/// Reusable polarized interferometer with an admitted, frequency-indexed arm cache.
///
/// Source grouping, native chromatic/patterned pupils and the analytic shutter
/// remain borrowed from the original exposure. Uniform arm elements are cached in
/// that owner's frequency-group order; observation performs no spectral search or
/// allocation. The immutable borrow prevents changing the pupil/source owner while
/// this receiver is live. A sample contains densities; a pixel collects W and J.
pub struct PreparedCoherentJonesInterferometer<'query, 'camera, 'source> {
    source: &'query PreparedCoherentJonesDiffractionExposure<'camera, 'source>,
    receiver: CoherentDiffractionInterferometer,
    arms: Vec<[CoherentJonesArm; 2]>,
}

impl<'camera, 'source> PreparedCoherentJonesDiffractionExposure<'camera, 'source> {
    /// Cache a complete passive arm spectrum on this exposure's sensor plane.
    ///
    /// Responses may arrive in any order, but exactly one must match each
    /// distinct authored frequency. Duplicate, missing or foreign frequencies
    /// refuse; equal-frequency source entries cannot acquire different optics.
    /// Relative common phases are retained, so sampled dispersive delays and
    /// filters change the original temporal interference. No passive loss is
    /// renormalized, and no unprovided frequency response is inferred.
    ///
    /// The byte limit covers the retained arm-array payload:
    /// frequency_count() * size_of::<[CoherentJonesArm; 2]>(), excluding the
    /// Vec header and allocator metadata. It is admitted before allocation.
    /// Source admission already bounds frequencies to 64 and each native arm to
    /// eight layers. Bounded O(N²) matching occurs only during preparation.
    /// Cancellation or failure publishes no receiver and never mutates the source.
    #[allow(clippy::float_cmp)] // Sensor plane and frequency identities are authored values.
    pub fn prepare_interferometer<'query>(
        &'query self,
        receiver: CoherentDiffractionInterferometer,
        responses: &[CoherentJonesInterferometerResponse],
        maximum_cached_arm_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        PreparedCoherentJonesInterferometer<'query, 'camera, 'source>,
        CoherentDiffractionError,
    > {
        poll(&mut cancelled)?;
        if receiver.sensor_distance_m() != self.optics.sensor_distance_m {
            return Err(CoherentDiffractionError::InvalidSensorPlane);
        }
        if responses.len() != self.frequency_count() {
            return Err(CoherentDiffractionError::InconsistentJonesArmSpectrum);
        }
        let quantity = "Jones interferometer cached arm bytes";
        let bytes = self
            .frequency_count()
            .checked_mul(core::mem::size_of::<[CoherentJonesArm; 2]>())
            .ok_or_else(|| limit(quantity, maximum_cached_arm_bytes))?;
        admit_count(quantity, bytes, maximum_cached_arm_bytes)?;
        for (index, response) in responses.iter().enumerate() {
            poll(&mut cancelled)?;
            if !response.frequency_hz.is_finite() || response.frequency_hz <= 0.0 {
                return Err(CoherentDiffractionError::InvalidFrequency);
            }
            if responses[..index]
                .iter()
                .any(|other| other.frequency_hz == response.frequency_hz)
                || !self
                    .optics
                    .groups
                    .iter()
                    .any(|group| group.frequency_hz == response.frequency_hz)
            {
                return Err(CoherentDiffractionError::InconsistentJonesArmSpectrum);
            }
        }
        poll(&mut cancelled)?;
        let mut arms = Vec::with_capacity(self.frequency_count());
        for group in &self.optics.groups {
            poll(&mut cancelled)?;
            let response = responses
                .iter()
                .find(|response| response.frequency_hz == group.frequency_hz)
                .ok_or(CoherentDiffractionError::InconsistentJonesArmSpectrum)?;
            arms.push(response.arms);
        }
        poll(&mut cancelled)?;
        Ok(PreparedCoherentJonesInterferometer {
            source: self,
            receiver,
            arms,
        })
    }
}

impl PreparedCoherentJonesInterferometer<'_, '_, '_> {
    /// Number of distinct original frequencies with an explicit arm response.
    #[must_use]
    pub fn frequency_count(&self) -> usize {
        self.arms.len()
    }

    /// Retained arm-array payload bytes, excluding Vec and allocator metadata.
    #[must_use]
    pub fn cached_arm_bytes(&self) -> usize {
        self.arms.len() * core::mem::size_of::<[CoherentJonesArm; 2]>()
    }

    /// Observe both full output coherencies on the original finite shutter.
    ///
    /// Source-column identity, exact active-arm work, DD destructive differences,
    /// native coherency admission and cancellation use the uniform receiver's
    /// shared implementation. The source's frequency-specific pupils precede
    /// these frequency-specific arm elements.
    pub fn observe(
        &self,
        limits: CoherentInterferometerExposureLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedCoherentInterferometerExposureSample, CoherentDiffractionError> {
        self.source.observe_interferometer_with_arms(
            self.receiver,
            InterferometerArms::Spectral(&self.arms),
            limits,
            cancelled,
        )
    }

    /// Collect the same chromatic arm fields over a finite relay pixel and shutter.
    ///
    /// All twelve moments share the original positive spatial quadrature and
    /// cross-frequency/cross-arm phase resolution. Visibility follows integrated
    /// complex mutual coherence. The native pixel budgets remain authoritative;
    /// no allocation or spectrum matching occurs during observation.
    pub fn integrate_pixel(
        &self,
        relay_offsets: CoherentDiffractionPixel,
        limits: CoherentDiffractionExposurePixelLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedCoherentInterferometerPixelExposure, CoherentDiffractionError> {
        self.source.integrate_interferometer_pixel_with_arms(
            self.receiver,
            InterferometerArms::Spectral(&self.arms),
            relay_offsets,
            limits,
            cancelled,
        )
    }
}
