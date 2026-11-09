//! Joint spatial, temporal and polarization integration of the two-arm field.

use super::super::super::super::{
    CoherentDiffractionError, CoherentDiffractionExposurePixelLimits,
    CoherentDiffractionInterferometer, CoherentDiffractionPixel,
    CoherentInterferometerExposureLimits, CoherentJonesArm, PolarizedDiffractionCoherency,
    PreparedCoherentJonesDiffractionExposure, interferometer::visibility, pixel::integrate_moments,
    poll, rounded,
};
use super::super::super::interferometer::energy;
use super::super::coherency::{energy as coherency_energy, physical};
use super::{InterferometerArms, InterferometerMoments};
use fs_math::{c64::C64, dd::Dd};

/// Two-port polarized energy collected over a finite relay pixel and shutter.
///
/// Coherency matrices and all four native Stokes components are spatially
/// integrated powers in watts, or shutter-integrated energies in joules. The
/// result can be depolarized by spatial averaging, temporal averaging or source
/// statistics; no single Jones field is assigned to that mixed observation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedCoherentInterferometerPixelExposure {
    /// Collected arm powers after the passive split and all authored arm optics.
    pub arm_power_w: [f64; 2],
    /// Complete average coherencies at the plus and minus output ports, in W.
    pub output_coherency_w: [PolarizedDiffractionCoherency; 2],
    /// Spatial integral of the shutter-averaged trace mutual coherence, in W.
    pub mutual_coherence_w: C64,
    /// Total-power contrast under an achromatic combiner-phase scan.
    pub fringe_visibility: f64,
    pub arm_energy_j: [f64; 2],
    /// Complete shutter-integrated output coherencies and signed Stokes, in J.
    pub output_coherency_j: [PolarizedDiffractionCoherency; 2],
    pub subdivisions: [usize; 2],
    pub quadrature_samples: usize,
    pub pupil_evaluations: usize,
    pub pair_contractions: usize,
}

impl PreparedCoherentJonesDiffractionExposure<'_, '_> {
    /// Collect polarized spectral interference over a physical detector pixel.
    ///
    /// The common tangent-coordinate `relay_offsets` translate both receiver
    /// coordinates through equal-area, unit-magnification relays. The original
    /// paraxial admission includes every translated corner of both arms. A
    /// physical lateral offset x on the sensor at axial distance d maps to x/d;
    /// the irradiance and area Jacobians cancel in the collected watts.
    ///
    /// At every positive spatial quadrature sample, each independent source
    /// bank propagates its original chromatic Jones fields through the passive
    /// arm elements. Fields combine before the analytic shutter contracts their
    /// full coherencies. All twelve DD moments share the same spatial rule;
    /// neither rounded point powers nor an averaged field replace this operation.
    /// Visibility is computed after integrating complex mutual coherence.
    ///
    /// Resolution includes the native pupil bandwidth and every retained
    /// cross-frequency/cross-arm sensor-phase gradient. Subdivision controls
    /// spatial resolution, not a certified error. The pupil, sample and temporal
    /// contraction counts are admitted before field evaluation. Fixed scratch
    /// allocates nothing; cancellation returns no partial observation and leaves
    /// preparation reusable. Uniform arm matrices are frequency-independent;
    /// authored chromatic pupil responses retain their original frequency scope.
    pub fn integrate_interferometer_pixel(
        &self,
        receiver: CoherentDiffractionInterferometer,
        arms: [CoherentJonesArm; 2],
        relay_offsets: CoherentDiffractionPixel,
        limits: CoherentDiffractionExposurePixelLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedCoherentInterferometerPixelExposure, CoherentDiffractionError> {
        self.integrate_interferometer_pixel_with_arms(
            receiver,
            InterferometerArms::Uniform(&arms),
            relay_offsets,
            limits,
            cancelled,
        )
    }

    pub(super) fn integrate_interferometer_pixel_with_arms(
        &self,
        receiver: CoherentDiffractionInterferometer,
        arms: InterferometerArms<'_>,
        relay_offsets: CoherentDiffractionPixel,
        limits: CoherentDiffractionExposurePixelLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedCoherentInterferometerPixelExposure, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        let (pupils_per_sample, pairs_per_sample) = self.jones_interferometer_work(
            receiver,
            CoherentInterferometerExposureLimits {
                maximum_pupil_evaluations: limits.spatial.maximum_pupil_evaluations,
                maximum_pair_contractions: limits.maximum_pair_contractions,
            },
        )?;
        let (subdivisions, quadrature_samples, pupil_evaluations, pair_contractions) =
            self.optics.interferometer_pixel_work(
                receiver,
                relay_offsets,
                limits,
                pupils_per_sample,
                pairs_per_sample,
            )?;
        let diagonal_indices = [0, 1, 2, 3, 6, 7];
        let mut active = [false; 6];
        let integrated: [Dd; 12] = integrate_moments(
            relay_offsets,
            subdivisions,
            &mut cancelled,
            |offset, cancelled| {
                let moments =
                    self.jones_interferometer_moments(receiver, arms, offset, cancelled)?;
                let values = flatten(moments);
                for (flag, index) in active.iter_mut().zip(diagonal_indices) {
                    *flag |= values[index] != Dd::ZERO;
                }
                Ok(values)
            },
        )?;
        // Positive spatial weights cannot make an observed population vanish.
        // Refuse a wholly underflowed collected power instead of declaring an
        // exact optical extinction. Genuine dark ports never set their flag.
        for (was_active, index) in active.into_iter().zip(diagonal_indices) {
            if was_active && integrated[index] == Dd::ZERO {
                return Err(CoherentDiffractionError::Unrepresentable);
            }
        }
        let powers = [integrated[0], integrated[1]];
        let outputs: [[Dd; 4]; 2] =
            std::array::from_fn(|port| std::array::from_fn(|axis| integrated[2 + 4 * port + axis]));
        let mutual = [integrated[10], integrated[11]];
        let duration = Dd::from_f64(self.shutter().effective_duration_s());
        let result = PolarizedCoherentInterferometerPixelExposure {
            arm_power_w: [rounded(powers[0])?, rounded(powers[1])?],
            output_coherency_w: [physical(outputs[0])?, physical(outputs[1])?],
            mutual_coherence_w: C64::new(rounded(mutual[0])?, rounded(mutual[1])?),
            fringe_visibility: visibility(powers, mutual)?,
            arm_energy_j: [energy(powers[0], duration)?, energy(powers[1], duration)?],
            output_coherency_j: [
                coherency_energy(outputs[0], duration)?,
                coherency_energy(outputs[1], duration)?,
            ],
            subdivisions,
            quadrature_samples,
            pupil_evaluations,
            pair_contractions,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
}

fn flatten(moments: InterferometerMoments) -> [Dd; 12] {
    std::array::from_fn(|i| match i {
        0..=1 => moments.arms[i],
        2..=5 => moments.outputs[0][i - 2],
        6..=9 => moments.outputs[1][i - 6],
        _ => moments.mutual[i - 10],
    })
}
