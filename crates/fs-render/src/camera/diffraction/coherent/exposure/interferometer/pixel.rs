//! Positive spatial integration of the original two-arm spectral receiver.
//!
//! A relay offset u translates both arm coordinates. The sensor-phase gradient
//! of arm a and optical frequency f is (2*pi/c)*d*f*(center_a+u). At either
//! endpoint of a pixel axis, the largest pair difference is the range of these
//! gradients over active arms and the two extreme optical frequencies. That
//! range is convex in u, so the endpoint maximum bounds every frequency pair
//! throughout the pixel, without an O(N²) bandwidth scan. Subtracting one common
//! reference gradient before finding the range preserves close-tone differences.

use super::super::super::interferometer::CoherentInterferometerPixelPower;
use super::super::{
    C_M_S, CoherentDiffractionExposurePixelLimits, CoherentDiffractionPixel, TAU, admit_coordinate,
    pixel,
};
use super::{
    CoherentDiffractionError, CoherentDiffractionInterferometer, Dd,
    PreparedCoherentDiffractionExposure, admit_count, energy, limit, poll, rounded, visibility,
};
use fs_math::c64::C64;

/// Spatially collected average power and finite-shutter energy at both ports.
/// The average-power coherences and visibility use the complete pixel integral.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentInterferometerPixelExposure {
    /// Collected watts, complex mutual coherence and exact spatial/pupil work.
    pub average_power: CoherentInterferometerPixelPower,
    /// Shutter-integrated arm energies after the passive input split, in joules.
    pub arm_energy_j: [f64; 2],
    /// Shutter-integrated plus/minus output energies, in joules.
    pub output_energy_j: [f64; 2],
    /// Total real frequency contractions over every quadrature sample and bank.
    pub pair_contractions: usize,
}

impl PreparedCoherentDiffractionExposure<'_, '_> {
    /// Collect a spectral interferogram over a finite relay pixel and shutter.
    ///
    /// `relay_offsets` translates both receiver coordinates by the same tangent
    /// offset. The complete translated rectangles, including inactive arms,
    /// must lie in the original paraxial domain. Unit-magnification, equal-area
    /// relays and the original nondispersive arm OPD convention remain.
    ///
    /// Arm and output powers and complex mutual coherence are integrated as DD
    /// moments through the native positive spatial rule. Each source bank's
    /// fields recombine before the analytic temporal contraction. Visibility is
    /// computed from the integrated coherence, so both spatial and temporal
    /// fringe washout are retained; it is not an average of point visibility.
    ///
    /// Resolution covers pupil bandwidth and every active cross-arm/frequency
    /// sensor chirp. It is a phase-span setting, not an error certificate. The
    /// exact sample, pupil and contraction counts are admitted before evaluating
    /// a field. Queries allocate nothing and cancellation leaves preparation
    /// reusable. Refinement changes spatial work, never introduces carrier or
    /// shutter time stepping.
    pub fn integrate_interferometer_pixel(
        &self,
        receiver: CoherentDiffractionInterferometer,
        relay_offsets: CoherentDiffractionPixel,
        limits: CoherentDiffractionExposurePixelLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentInterferometerPixelExposure, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        let per_sample =
            receiver.pupil_evaluations(self.waves, limits.spatial.maximum_pupil_evaluations)?;
        let active = receiver
            .arm_amplitudes()
            .iter()
            .filter(|&&a| a != Dd::ZERO)
            .count();
        let per_sample_pairs =
            self.interferometer_work(active, limits.maximum_pair_contractions)?;
        let (subdivisions, quadrature_samples, pupil_evaluations, pair_contractions) = self
            .interferometer_pixel_work(
                receiver,
                relay_offsets,
                limits,
                per_sample,
                per_sample_pairs,
            )?;
        let mut nonzero = [false; 4];
        let moments = pixel::integrate_moments(
            relay_offsets,
            subdivisions,
            &mut cancelled,
            |offset, cancelled| {
                let value = self.interferometer_moments(receiver, offset, cancelled)?;
                for (flag, &power) in nonzero.iter_mut().zip(&value) {
                    *flag |= power != Dd::ZERO;
                }
                Ok(value)
            },
        )?;
        for (active, &power) in nonzero.iter().zip(&moments) {
            if *active && power == Dd::ZERO {
                return Err(CoherentDiffractionError::Unrepresentable);
            }
        }
        let arms = [moments[0], moments[1]];
        let outputs = [moments[2], moments[3]];
        let mutual = [moments[4], moments[5]];
        let duration = Dd::from_f64(self.shutter.duration_s());
        let result = CoherentInterferometerPixelExposure {
            average_power: CoherentInterferometerPixelPower {
                arm_power_w: [rounded(arms[0])?, rounded(arms[1])?],
                output_power_w: [rounded(outputs[0])?, rounded(outputs[1])?],
                mutual_coherence_w: C64::new(rounded(mutual[0])?, rounded(mutual[1])?),
                fringe_visibility: visibility(arms, mutual)?,
                subdivisions,
                quadrature_samples,
                pupil_evaluations,
            },
            arm_energy_j: [energy(arms[0], duration)?, energy(arms[1], duration)?],
            output_energy_j: [energy(outputs[0], duration)?, energy(outputs[1], duration)?],
            pair_contractions,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    // Shared by scalar and Jones interferograms: the source/contract counts
    // differ, while the original pupil, arm geometry and sensor phases do not.
    #[allow(clippy::float_cmp)] // The prepared sensor plane is an authored value.
    pub(in super::super) fn interferometer_pixel_work(
        &self,
        receiver: CoherentDiffractionInterferometer,
        relay_offsets: CoherentDiffractionPixel,
        limits: CoherentDiffractionExposurePixelLimits,
        pupil_evaluations_per_sample: usize,
        pair_contractions_per_sample: usize,
    ) -> Result<([usize; 2], usize, usize, usize), CoherentDiffractionError> {
        if receiver.sensor_distance_m() != self.sensor_distance_m {
            return Err(CoherentDiffractionError::InvalidSensorPlane);
        }
        let centers = receiver.coordinates();
        let lower = relay_offsets.minimum_tangent();
        let upper = relay_offsets.maximum_tangent();
        for x in [lower[0], upper[0]] {
            for y in [lower[1], upper[1]] {
                for center in centers {
                    admit_coordinate([
                        (Dd::from_f64(center[0]) + Dd::from_f64(x)).to_f64(),
                        (Dd::from_f64(center[1]) + Dd::from_f64(y)).to_f64(),
                    ])?;
                }
            }
        }
        let mut maximum_step = [1.0; 2];
        if pupil_evaluations_per_sample != 0 {
            let first = self
                .groups
                .first()
                .ok_or(CoherentDiffractionError::Unrepresentable)?;
            let mut minimum_frequency = first.optical_frequency_hz;
            let mut maximum_frequency = minimum_frequency;
            for group in &self.groups[1..] {
                if group.optical_frequency_hz.lt(minimum_frequency) {
                    minimum_frequency = group.optical_frequency_hz;
                }
                if maximum_frequency.lt(group.optical_frequency_hz) {
                    maximum_frequency = group.optical_frequency_hz;
                }
            }
            let frequency_span = maximum_frequency - minimum_frequency;
            let amplitudes = receiver.arm_amplitudes();
            let reference = usize::from(amplitudes[0] == Dd::ZERO);
            let pupil_bandwidth = Dd::from_f64(2.0) / Dd::from_f64(self.minimum_scale);
            for axis in 0..2 {
                let mut maximum_range = Dd::ZERO;
                for offset in [lower[axis], upper[axis]] {
                    let mut low = Dd::ZERO;
                    let mut high = Dd::ZERO;
                    for arm in 0..2 {
                        if amplitudes[arm] == Dd::ZERO {
                            continue;
                        }
                        let coordinate = Dd::from_f64(centers[arm][axis]);
                        let base = minimum_frequency
                            * (coordinate - Dd::from_f64(centers[reference][axis]));
                        for value in [
                            base,
                            base + frequency_span * (coordinate + Dd::from_f64(offset)),
                        ] {
                            if value.lt(low) {
                                low = value;
                            }
                            if high.lt(value) {
                                high = value;
                            }
                        }
                    }
                    let range = high - low;
                    if maximum_range.lt(range) {
                        maximum_range = range;
                    }
                }
                let chirp_bandwidth = (maximum_range / Dd::from_f64(C_M_S))
                    * TAU
                    * Dd::from_f64(self.sensor_distance_m);
                maximum_step[axis] = limits.spatial.maximum_phase_span_rad
                    / rounded(pupil_bandwidth + chirp_bandwidth)?;
            }
        }
        let (subdivisions, samples, pupil_evaluations) = pixel::pixel_work(
            relay_offsets,
            limits.spatial,
            maximum_step,
            pupil_evaluations_per_sample,
        )?;
        let pair_contractions = samples
            .checked_mul(pair_contractions_per_sample)
            .ok_or_else(|| {
                limit(
                    "pixel interferometer pair contractions",
                    limits.maximum_pair_contractions,
                )
            })?;
        admit_count(
            "pixel interferometer pair contractions",
            pair_contractions,
            limits.maximum_pair_contractions,
        )?;
        Ok((subdivisions, samples, pupil_evaluations, pair_contractions))
    }
}
