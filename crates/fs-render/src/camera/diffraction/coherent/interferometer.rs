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
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentInterferometerSample, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        admit_sensor(self.camera, receiver.sensor_distance_m)?;
        let active = receiver
            .arm_power_fractions
            .iter()
            .filter(|&&p| p != 0.0)
            .count();
        let evaluations =
            self.waves
                .checked_mul(active)
                .ok_or(CoherentDiffractionError::Limit {
                    quantity: "interferometer pupil evaluations",
                    requested: usize::MAX,
                    maximum: maximum_pupil_evaluations,
                })?;
        admit_count(
            "interferometer pupil evaluations",
            evaluations,
            maximum_pupil_evaluations,
        )?;
        let amplitudes = receiver.arm_power_fractions.map(|p| Dd::from_f64(p).sqrt());
        // Query the original medium wavelength, not a group-delay index.
        let medium_wavelength_nm = self.camera.pupil_wavelength_nm(self.wavelength_nm)?;
        let radii = receiver.coordinates.map(|at| {
            Dd::from_f64(at[0]) * Dd::from_f64(at[0]) + Dd::from_f64(at[1]) * Dd::from_f64(at[1])
        });
        // Keep independently authored phases separate until each has been
        // reduced; a large OPD must not erase a small physical phase offset.
        let phases = [
            length_phase(
                Dd::from_f64(receiver.sensor_distance_m)
                    * (radii[1] - radii[0])
                    * Dd::from_f64(0.5),
                medium_wavelength_nm,
            )?,
            length_phase(
                Dd::from_f64(receiver.optical_path_difference_m),
                self.wavelength_nm,
            )?,
            Dd::from_f64(receiver.relative_phase_rad),
        ];
        let mut rotation = [Dd::ONE, Dd::ZERO];
        for phase in phases {
            rounded(phase)?;
            let (sine, cosine) = det::sin_cos_dd(phase);
            rotation = mul(rotation, [cosine, sine]);
        }
        let mut arms = [Dd::ZERO; 2];
        let mut outputs = [Dd::ZERO; 2];
        let mut mutual = [Dd::ZERO; 2];
        for bank in core::iter::once(self.illumination.mean)
            .chain(self.illumination.factors.iter().map(|factor| factor.waves))
        {
            poll(&mut cancelled)?;
            let mut fields = [[Dd::ZERO; 2]; 2];
            for arm in 0..2 {
                if receiver.arm_power_fractions[arm] != 0.0 {
                    fields[arm] = self
                        .field(bank, receiver.coordinates[arm], &mut cancelled)?
                        .map(|value| value * amplitudes[arm]);
                }
            }
            fields[1] = mul(fields[1], rotation);
            for arm in 0..2 {
                arms[arm] = arms[arm] + square(fields[arm])?;
            }
            let cross = mul(fields[0], [fields[1][0], -fields[1][1]]);
            for axis in 0..2 {
                mutual[axis] = mutual[axis] + cross[axis];
            }
            for (port, sign) in [Dd::ONE, -Dd::ONE].into_iter().enumerate() {
                let field = std::array::from_fn(|axis| fields[0][axis] + sign * fields[1][axis]);
                outputs[port] = outputs[port] + square(field)? * Dd::from_f64(0.5);
            }
        }
        let total = arms[0] + arms[1];
        let visibility = if total == Dd::ZERO {
            0.0
        } else {
            // Normalize before squaring to avoid overflow in bright fields.
            let re = mutual[0] / total;
            let im = mutual[1] / total;
            rounded((re * re + im * im).sqrt() * Dd::from_f64(2.0))?
        };
        if visibility > 1.0 + 16.0 * f64::EPSILON {
            return Err(CoherentDiffractionError::Unrepresentable);
        }
        let result = CoherentInterferometerSample {
            arm_power_density_w: [rounded(arms[0])?, rounded(arms[1])?],
            output_power_density_w: [rounded(outputs[0])?, rounded(outputs[1])?],
            mutual_coherence_density_w: C64::new(rounded(mutual[0])?, rounded(mutual[1])?),
            fringe_visibility: visibility.min(1.0),
            pupil_evaluations: evaluations,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
}

// Reduce the optical length before the irrational radian conversion. Using
// lambda_nm * f64(1e-9) first would inject the rounded unit-conversion error
// into every cycle, producing a false dark-port signal even for a 1 m OPD.
// The exact integer metre-to-nanometre factor preserves the authored units;
// reducing both retained components before division also keeps large paths
// from accumulating the wavelength-division or 2*pi approximation errors.
fn length_phase(length_m: Dd, wavelength_nm: f64) -> Result<Dd, CoherentDiffractionError> {
    let length_nm = length_m * Dd::from_f64(1e9);
    rounded(length_nm)?;
    let remainder_nm = Dd::from_pair(length_nm.hi % wavelength_nm, length_nm.lo % wavelength_nm);
    let tau = Dd::from_pair(core::f64::consts::TAU, 2.449_293_598_294_706_4e-16);
    let phase = (remainder_nm / Dd::from_f64(wavelength_nm)) * tau;
    rounded(phase)?;
    Ok(phase)
}

fn mul(a: [Dd; 2], b: [Dd; 2]) -> [Dd; 2] {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
