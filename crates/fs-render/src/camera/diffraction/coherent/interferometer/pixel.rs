//! Finite equal-area relay pixels of the original two-arm receiver.

use super::super::{
    CoherentDiffractionError, CoherentDiffractionPixel, CoherentDiffractionPixelLimits, Dd,
    PreparedCoherentDiffraction, admit_coordinate,
    pixel::{integrate_moments, pixel_work},
    poll, rounded,
};
use super::{
    CoherentDiffractionInterferometer, admit_sensor, length_phase, mul, rotation, visibility,
};
use fs_math::c64::C64;

/// Spatially integrated monochromatic two-arm observation, in watts.
/// The footprint translates both sensor coordinates by the same relay offset.
/// Unit magnification and equal tangent-area measure are required in both arms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentInterferometerPixelPower {
    /// Collected arm powers after passive pickoff and before recombination.
    pub arm_power_w: [f64; 2],
    /// Collected plus/minus port powers, formed from their complete fields.
    pub output_power_w: [f64; 2],
    /// Integral of W01 = <E0 conj(E1)> after propagation/interface phase.
    pub mutual_coherence_w: C64,
    /// Integrated fringe contrast 2*|integral W01|/(P0+P1), or zero in darkness.
    /// Spatially varying fringes can reduce this below every point visibility.
    pub fringe_visibility: f64,
    /// Subcells on the original pupil axes, or zero when no arm observes waves.
    pub subdivisions: [usize; 2],
    /// Number of positive area-quadrature samples used.
    pub quadrature_samples: usize,
    /// Original pupil-transfer evaluations across every active arm and sample.
    pub pupil_evaluations: usize,
}

impl PreparedCoherentDiffraction<'_, '_> {
    /// Collect both interferometer outputs over a common finite relay pixel.
    ///
    /// `relay_offsets` are tangent OFFSETS from each receiver coordinate, not
    /// absolute sensor positions. The complete translated rectangles must stay
    /// inside the original paraxial domain, including inactive arms. Equal-area,
    /// unit-magnification relays bring corresponding samples to the combiner.
    /// For a physical sensor at distance d, lateral offsets are divided by d;
    /// the irradiance and area Jacobians cancel in the collected power.
    ///
    /// Each mean/covariance bank interferes before positive quadrature weights
    /// multiply its port and arm powers. Mutual coherence retains its complex
    /// sign through the same DD integration. No rounded point-power difference
    /// is used to obtain a dark port, and independent factors never interfere.
    ///
    /// Resolution includes both pupil intensity bandwidth and the relative
    /// sensor chirp: 2/scale + |k_medium*d*(center1-center0)| per axis. The
    /// complete sample and active-arm pupil work is admitted before evaluation.
    /// Refining the phase span assesses quadrature convergence; it is not an
    /// optical accuracy certificate. Evaluation allocates nothing, polls each
    /// sample/wave, and leaves the preparation reusable after cancellation.
    pub fn integrate_interferometer_pixel(
        &self,
        receiver: CoherentDiffractionInterferometer,
        relay_offsets: CoherentDiffractionPixel,
        limits: CoherentDiffractionPixelLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentInterferometerPixelPower, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        admit_sensor(self.camera, receiver.sensor_distance_m())?;
        let centers = receiver.coordinates();
        let lower = relay_offsets.minimum_tangent();
        let upper = relay_offsets.maximum_tangent();
        for x in [lower[0], upper[0]] {
            for y in [lower[1], upper[1]] {
                for coordinate in translated(centers, [x, y]) {
                    admit_coordinate(coordinate)?;
                }
            }
        }
        let amplitudes = receiver.arm_amplitudes();
        let per_sample =
            receiver.pupil_evaluations(self.waves, limits.maximum_pupil_evaluations)?;
        let interfering = per_sample != 0 && amplitudes.iter().all(|&a| a != Dd::ZERO);
        let medium_wavelength_nm = self.camera.pupil_wavelength_nm(self.wavelength_nm)?;
        // Common offsets cancel the quadratic term exactly. Keeping this
        // linear path difference apart from the centre phase avoids subtracting
        // independently rounded chirps at every quadrature sample.
        let path_gradient: [Dd; 2] = if interfering {
            std::array::from_fn(|axis| {
                Dd::from_f64(receiver.sensor_distance_m())
                    * (Dd::from_f64(centers[1][axis]) - Dd::from_f64(centers[0][axis]))
            })
        } else {
            [Dd::ZERO; 2]
        };
        let mut maximum_step = [1.0; 2];
        if per_sample != 0 {
            let pupil_bandwidth = Dd::from_f64(2.0) / Dd::from_f64(self.pupil.scale);
            let tau = Dd::from_pair(core::f64::consts::TAU, 2.449_293_598_294_706_4e-16);
            for axis in 0..2 {
                let chirp_bandwidth = path_gradient[axis].abs() * Dd::from_f64(1e9)
                    / Dd::from_f64(medium_wavelength_nm)
                    * tau;
                maximum_step[axis] =
                    limits.maximum_phase_span_rad / rounded(pupil_bandwidth + chirp_bandwidth)?;
            }
        }
        let (subdivisions, quadrature_samples, pupil_evaluations) =
            pixel_work(relay_offsets, limits, maximum_step, per_sample)?;
        let center_rotation = if interfering {
            receiver.relative_rotation(self.wavelength_nm, medium_wavelength_nm, centers)?
        } else {
            [Dd::ONE, Dd::ZERO]
        };
        let moments = integrate_moments(
            relay_offsets,
            subdivisions,
            &mut cancelled,
            |offset, cancelled| {
                let arm_rotation = if interfering {
                    let path = path_gradient[0] * Dd::from_f64(offset[0])
                        + path_gradient[1] * Dd::from_f64(offset[1]);
                    mul(
                        center_rotation,
                        rotation([length_phase(path, medium_wavelength_nm)?])?,
                    )
                } else {
                    center_rotation
                };
                self.interferometer_moments(
                    translated(centers, offset),
                    amplitudes,
                    arm_rotation,
                    cancelled,
                )
            },
        )?;
        let result = CoherentInterferometerPixelPower {
            arm_power_w: [rounded(moments[0])?, rounded(moments[1])?],
            output_power_w: [rounded(moments[2])?, rounded(moments[3])?],
            mutual_coherence_w: C64::new(rounded(moments[4])?, rounded(moments[5])?),
            fringe_visibility: visibility([moments[0], moments[1]], [moments[4], moments[5]])?,
            subdivisions,
            quadrature_samples,
            pupil_evaluations,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
}

fn translated(centers: [[f64; 2]; 2], offset: [f64; 2]) -> [[f64; 2]; 2] {
    centers.map(|center| {
        std::array::from_fn(|axis| {
            (Dd::from_f64(center[axis]) + Dd::from_f64(offset[axis])).to_f64()
        })
    })
}
