//! Two-arm detection of the original camera's mutual spatial coherence.
//!
//! Fields, including each shared covariance column, pass through a passive
//! pickoff and a balanced lossless combiner. The outputs are (a+b)/sqrt(2)
//! and (a-b)/sqrt(2); their powers sum to the admitted two-arm input power.
//! A fixed output phase convention absorbs the combiner's quadrature into
//! the caller's explicit relative phase. No coherence is inferred from two
//! intensity samples. This is a scalar, monochromatic, paraxial receiver,
//! not propagation through a resolved interferometer or detector electronics.

use super::{
    CoherentDiffractionError, Dd, PreparedCoherentDiffraction, admit_coordinate, admit_count,
    exposure::admit_sensor, poll, rounded, square,
};
use fs_math::{c64::C64, det};

mod pixel;
pub use pixel::CoherentInterferometerPixelPower;
mod linewidth;
pub use linewidth::InterferometerTemporalCoherence;
pub(super) use linewidth::TemporalWeights;

/// Passive sampling of two positions on one physical sensor plane.
///
/// The arm powers are fractions of the original clear-pupil source power.
/// Their sum may not exceed one: coincident samples therefore cannot duplicate
/// source power. Separate source banks remain uncorrelated through the receiver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentDiffractionInterferometer {
    coordinates: [[f64; 2]; 2],
    sensor_distance_m: f64,
    arm_power_fractions: [f64; 2],
    optical_path_difference_m: f64,
    relative_phase_rad: f64,
}

impl CoherentDiffractionInterferometer {
    /// Positive path difference delays arm 1 relative to arm 0. With the
    /// original exp(-i omega t) convention its phase is +2*pi*OPD/lambda_vac.
    /// OPD is optical length (including any arm medium's phase index), not
    /// geometric length. `relative_phase_rad` carries other authored interface
    /// or combiner phases. Neither quantity is inferred from the scene.
    pub fn try_new(
        coordinates: [[f64; 2]; 2],
        sensor_distance_m: f64,
        arm_power_fractions: [f64; 2],
        optical_path_difference_m: f64,
        relative_phase_rad: f64,
    ) -> Result<Self, CoherentDiffractionError> {
        for coordinate in coordinates {
            admit_coordinate(coordinate)?;
        }
        if !sensor_distance_m.is_finite() || sensor_distance_m <= 0.0 {
            return Err(CoherentDiffractionError::InvalidSensorPlane);
        }
        if arm_power_fractions
            .iter()
            .any(|p| !p.is_finite() || !(0.0..=1.0).contains(p))
            || Dd::ONE
                .lt(Dd::from_f64(arm_power_fractions[0]) + Dd::from_f64(arm_power_fractions[1]))
            || !optical_path_difference_m.is_finite()
            || !relative_phase_rad.is_finite()
        {
            return Err(CoherentDiffractionError::InvalidInterferometer);
        }
        Ok(Self {
            coordinates,
            sensor_distance_m,
            arm_power_fractions,
            optical_path_difference_m,
            relative_phase_rad,
        })
    }

    pub(super) const fn coordinates(self) -> [[f64; 2]; 2] {
        self.coordinates
    }

    pub(super) const fn sensor_distance_m(self) -> f64 {
        self.sensor_distance_m
    }

    pub(super) const fn optical_path_difference_m(self) -> f64 {
        self.optical_path_difference_m
    }

    pub(super) const fn relative_phase_rad(self) -> f64 {
        self.relative_phase_rad
    }

    pub(super) fn arm_amplitudes(self) -> [Dd; 2] {
        self.arm_power_fractions
            .map(|power| Dd::from_f64(power).sqrt())
    }

    pub(super) fn pupil_evaluations(
        self,
        waves: usize,
        maximum: usize,
    ) -> Result<usize, CoherentDiffractionError> {
        let active = self
            .arm_power_fractions
            .iter()
            .filter(|&&p| p != 0.0)
            .count();
        let evaluations = waves
            .checked_mul(active)
            .ok_or(CoherentDiffractionError::Limit {
                quantity: "interferometer pupil evaluations",
                requested: usize::MAX,
                maximum,
            })?;
        admit_count("interferometer pupil evaluations", evaluations, maximum)?;
        Ok(evaluations)
    }

    // Monochromatic pupil fields omit the common sensor chirp. Restore its
    // relative value at the two actual positions before recombination.
    pub(super) fn relative_rotation(
        self,
        wavelength_nm: f64,
        medium_wavelength_nm: f64,
        coordinates: [[f64; 2]; 2],
    ) -> Result<[Dd; 2], CoherentDiffractionError> {
        let radii = coordinates.map(|at| {
            Dd::from_f64(at[0]) * Dd::from_f64(at[0]) + Dd::from_f64(at[1]) * Dd::from_f64(at[1])
        });
        rotation([
            length_phase(
                Dd::from_f64(self.sensor_distance_m) * (radii[1] - radii[0]) * Dd::from_f64(0.5),
                medium_wavelength_nm,
            )?,
            length_phase(Dd::from_f64(self.optical_path_difference_m), wavelength_nm)?,
            Dd::from_f64(self.relative_phase_rad),
        ])
    }
}

/// Ensemble powers per tangent-coordinate area at a two-arm receiver.
///
/// The two input coordinates represent corresponding equal-area samples after
/// unit-magnification arm relays. A finite aperture or detector pixel requires
/// spatial integration; these point densities are not collected watts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentInterferometerSample {
    /// Arm powers after the passive input split and before recombination.
    pub arm_power_density_w: [f64; 2],
    /// Positive outputs in the declared plus/minus port convention.
    pub output_power_density_w: [f64; 2],
    /// W01 = <E0 conj(E1)> after the relative propagation/interface phase.
    pub mutual_coherence_density_w: C64,
    /// Fringe contrast 2*|W01|/(P0+P1) for scanning the relative phase.
    /// Zero when both arms are dark. Unequal powers reduce visibility.
    pub fringe_visibility: f64,
    /// Original pupil-transfer evaluations actually performed.
    pub pupil_evaluations: usize,
}

impl PreparedCoherentDiffraction<'_, '_> {
    /// Recombine the original phase-aware fields at two sensor positions.
    ///
    /// The common Fresnel sensor chirp omitted from monochromatic intensity
    /// queries is restored as a RELATIVE phase before the arms interfere.
    /// Independent covariance columns add powers only after each column's two
    /// fields combine. Dark-port power is squared from the retained field
    /// difference, never obtained by subtracting rounded intensities.
    ///
    /// Work is O(waves + factors), bounded before any pupil evaluations, with
    /// constant scratch and no allocation. Cancelled/refused observations do
    /// not alter the reusable preparation.
    pub fn observe_interferometer(
        &self,
        receiver: CoherentDiffractionInterferometer,
        maximum_pupil_evaluations: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentInterferometerSample, CoherentDiffractionError> {
        self.observe_interferometer_with_coherence(
            receiver, None, maximum_pupil_evaluations, cancelled,
        )
    }

    /// Stationary finite-linewidth port powers using the original spatial fields.
    ///
    /// The explicit narrowband source multiplies cross-arm coherence by its
    /// analytic temporal envelope at the authored group delay. Each bank's
    /// output is a positive mixture of its coherent port and independent-arm
    /// powers. Passive arm energy is unchanged, including at a dark fringe.
    /// This ensemble mean does not model finite-shutter noise or infer delays
    /// from the phase OPD. See [`InterferometerTemporalCoherence`].
    pub fn observe_interferometer_with_linewidth(
        &self,
        receiver: CoherentDiffractionInterferometer,
        coherence: InterferometerTemporalCoherence,
        maximum_pupil_evaluations: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentInterferometerSample, CoherentDiffractionError> {
        self.observe_interferometer_with_coherence(
            receiver, Some(coherence), maximum_pupil_evaluations, cancelled,
        )
    }

    fn observe_interferometer_with_coherence(
        &self,
        receiver: CoherentDiffractionInterferometer,
        coherence: Option<InterferometerTemporalCoherence>,
        maximum_pupil_evaluations: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentInterferometerSample, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        let temporal = if let Some(coherence) = coherence {
            coherence.admit(self.wavelength_nm)?;
            coherence.weights([0.0; 2])?
        } else {
            TemporalWeights::COHERENT
        };
        admit_sensor(self.camera, receiver.sensor_distance_m)?;
        let evaluations = receiver.pupil_evaluations(self.waves, maximum_pupil_evaluations)?;
        let amplitudes = receiver.arm_amplitudes();
        // Query the original medium wavelength, not a group-delay index.
        let medium_wavelength_nm = self.camera.pupil_wavelength_nm(self.wavelength_nm)?;
        let rotation = receiver.relative_rotation(
            self.wavelength_nm,
            medium_wavelength_nm,
            receiver.coordinates,
        )?;
        let moments = self.interferometer_moments(
            receiver.coordinates,
            amplitudes,
            rotation,
            temporal,
            &mut cancelled,
        )?;
        let visibility = visibility([moments[0], moments[1]], [moments[4], moments[5]])?;
        let result = CoherentInterferometerSample {
            arm_power_density_w: [rounded(moments[0])?, rounded(moments[1])?],
            output_power_density_w: [rounded(moments[2])?, rounded(moments[3])?],
            mutual_coherence_density_w: C64::new(rounded(moments[4])?, rounded(moments[5])?),
            fringe_visibility: visibility,
            pupil_evaluations: evaluations,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    // Preserve these six DD moments through finite detector integration.
    // Rounding a point sample's powers first would lose weak output signals.
    fn interferometer_moments(
        &self,
        coordinates: [[f64; 2]; 2],
        amplitudes: [Dd; 2],
        rotation: [Dd; 2],
        temporal: TemporalWeights,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Dd; 6], CoherentDiffractionError> {
        let mut arms = [Dd::ZERO; 2];
        let mut outputs = [Dd::ZERO; 2];
        let mut mutual = [Dd::ZERO; 2];
        let combiner_amplitude = Dd::from_f64(0.5).sqrt();
        for bank in core::iter::once(self.illumination.mean)
            .chain(self.illumination.factors.iter().map(|factor| factor.waves))
        {
            poll(cancelled)?;
            let mut fields = [[Dd::ZERO; 2]; 2];
            for arm in 0..2 {
                if amplitudes[arm] != Dd::ZERO {
                    fields[arm] = self
                        .field(bank, coordinates[arm], cancelled)?
                        .map(|value| value * amplitudes[arm]);
                }
            }
            fields[1] = mul(fields[1], rotation);
            let powers = [square(fields[0])?, square(fields[1])?];
            for arm in 0..2 {
                arms[arm] = arms[arm] + powers[arm];
            }
            let independent = powers[0] * (temporal.incoherent * Dd::from_f64(0.5))
                + powers[1] * (temporal.incoherent * Dd::from_f64(0.5));
            let cross = mul(fields[0], [fields[1][0], -fields[1][1]]);
            for axis in 0..2 {
                mutual[axis] = mutual[axis] + cross[axis] * temporal.coherent;
            }
            for (port, sign) in [Dd::ONE, -Dd::ONE].into_iter().enumerate() {
                let field = std::array::from_fn(|axis| {
                    ((fields[0][axis] + sign * fields[1][axis]) * combiner_amplitude)
                        * temporal.field_scale
                });
                outputs[port] = outputs[port] + (square(field)? + independent);
            }
        }
        Ok([
            arms[0], arms[1], outputs[0], outputs[1], mutual[0], mutual[1],
        ])
    }
}

// Reduce the optical length before the irrational radian conversion. Using
// lambda_nm * f64(1e-9) first would inject the rounded unit-conversion error
// into every cycle, producing a false dark-port signal even for a 1 m OPD.
// The exact integer metre-to-nanometre factor preserves the authored units;
// reducing both retained components before division also keeps large paths
// from accumulating the wavelength-division or 2*pi approximation errors.
pub(super) fn length_phase(
    length_m: Dd,
    wavelength_nm: f64,
) -> Result<Dd, CoherentDiffractionError> {
    let length_nm = length_m * Dd::from_f64(1e9);
    rounded(length_nm)?;
    let remainder_nm = Dd::from_pair(length_nm.hi % wavelength_nm, length_nm.lo % wavelength_nm);
    let tau = Dd::from_pair(core::f64::consts::TAU, 2.449_293_598_294_706_4e-16);
    let phase = (remainder_nm / Dd::from_f64(wavelength_nm)) * tau;
    rounded(phase)?;
    Ok(phase)
}

pub(super) fn mul(a: [Dd; 2], b: [Dd; 2]) -> [Dd; 2] {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

fn rotation(phases: impl IntoIterator<Item = Dd>) -> Result<[Dd; 2], CoherentDiffractionError> {
    let mut rotation = [Dd::ONE, Dd::ZERO];
    // Independently authored phases remain separate until each reduction, so
    // a large optical path cannot swallow a tiny interface phase.
    for phase in phases {
        rounded(phase)?;
        let (sine, cosine) = det::sin_cos_dd(phase);
        rotation = mul(rotation, [cosine, sine]);
    }
    Ok(rotation)
}

pub(super) fn visibility(arms: [Dd; 2], mutual: [Dd; 2]) -> Result<f64, CoherentDiffractionError> {
    let scale = if arms[0].lt(arms[1]) {
        arms[1]
    } else {
        arms[0]
    };
    let visibility = if scale == Dd::ZERO {
        0.0
    } else {
        // Every returned arm and port can be finite even when the sum of
        // the two arm powers exceeds binary64. Normalize before adding too.
        let total = arms[0] / scale + arms[1] / scale;
        let re = (mutual[0] / scale) / total;
        let im = (mutual[1] / scale) / total;
        rounded((re * re + im * im).sqrt() * Dd::from_f64(2.0))?
    };
    if visibility > 1.0 + 16.0 * f64::EPSILON {
        return Err(CoherentDiffractionError::Unrepresentable);
    }
    Ok(visibility.min(1.0))
}
