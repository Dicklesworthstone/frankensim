//! Analytic finite-exposure detection after a passive two-arm combiner.
//!
//! Each original covariance column remains one shared source across frequencies
//! and arms. Fields combine before the temporal Gram contraction, preserving a
//! weak dark port. The phase uses authored Hz, the original camera phase index,
//! and optical lengths; a rounded wavelength is never inverted for propagation.

use super::super::interferometer::{CoherentDiffractionInterferometer, mul, visibility};
use super::{
    C_M_S, CoherentDiffractionError, Dd, MAX_COHERENT_DIFFRACTION_TONES,
    PreparedCoherentDiffractionExposure, TAU, admit_count, limit, poll, rounded,
};
use fs_math::{c64::C64, det};

mod pixel;
pub use pixel::CoherentInterferometerPixelExposure;

type Fields = [[Dd; 2]; MAX_COHERENT_DIFFRACTION_TONES];

/// Work limits for one allocation-free spectral receiver observation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoherentInterferometerExposureLimits {
    /// Original pupil-wave evaluations across active arms and all source banks.
    pub maximum_pupil_evaluations: usize,
    /// Prepared temporal work across every source bank. If S is one scalar
    /// power contraction's work and C one complex mutual contraction's work,
    /// two active scalar arms use B*(4*S+C); one uses B*S, for B source banks.
    /// The short-band path counts its retained mode projections and products;
    /// pairwise fallback has S=N*(N+1)/2 and C=N*(2*N-1), for N frequencies.
    /// Full Jones coherency uses the larger counts documented by that query.
    /// A closed receiver or a preparation with no waves uses zero contractions.
    pub maximum_pair_contractions: usize,
}

/// Shutter-averaged ensemble densities at a two-arm receiver.
///
/// Coordinates denote corresponding equal-area samples after unit-magnification
/// relays. These are densities per tangent-coordinate area, not collected watts
/// or joules. A finite detector requires spatial integration. Visibility is the
/// contrast under an achromatic combiner-phase scan, not a broadband OPD scan.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentInterferometerExposureSample {
    pub arm_power_density_w: [f64; 2],
    pub output_power_density_w: [f64; 2],
    /// Window average of E0(t)*conj(E1(t)), after the relative arm phase.
    pub mutual_coherence_density_w: C64,
    pub fringe_visibility: f64,
    pub arm_energy_density_j: [f64; 2],
    pub output_energy_density_j: [f64; 2],
    pub pupil_evaluations: usize,
    pub pair_contractions: usize,
}

impl PreparedCoherentDiffractionExposure<'_, '_> {
    /// Recombine the original spectral fields and integrate their exact beats.
    ///
    /// The receiver must name this preparation's sensor plane. Its OPD is one
    /// shared, nondispersive optical length; no extra arm medium or coherence
    /// law is inferred. Camera propagation uses each admitted phase index, while
    /// OPD uses the original vacuum frequency. Both are referenced to the same
    /// pupil source clock and the original shutter, including negative epochs.
    ///
    /// Identical frequencies coalesce before any products. Each source bank's
    /// arm fields combine before squaring, rather than subtracting bright powers
    /// to obtain a dark port. Powers and mutual coherence use the same prepared
    /// short-band modes when available, preserving their common Gram structure.
    /// Fixed scratch and preflighted temporal/source work require no carrier
    /// or shutter time stepping. The pairwise fallback costs O(B*N^2 + waves). This is
    /// a scalar paraxial, cycle-averaged detector, not subcycle electric-field
    /// sampling, a polarization model, detector noise, or an error certificate.
    #[allow(clippy::float_cmp)] // The prepared sensor plane is an authored value.
    pub fn observe_interferometer(
        &self,
        receiver: CoherentDiffractionInterferometer,
        limits: CoherentInterferometerExposureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentInterferometerExposureSample, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        if receiver.sensor_distance_m() != self.sensor_distance_m {
            return Err(CoherentDiffractionError::InvalidSensorPlane);
        }
        let evaluations =
            receiver.pupil_evaluations(self.waves, limits.maximum_pupil_evaluations)?;
        let amplitudes = receiver.arm_amplitudes();
        let active = amplitudes.iter().filter(|&&a| a != Dd::ZERO).count();
        let contractions = self.interferometer_work(active, limits.maximum_pair_contractions)?;
        let moments = self.interferometer_moments(receiver, [0.0; 2], &mut cancelled)?;
        let arms = [moments[0], moments[1]];
        let outputs = [moments[2], moments[3]];
        let mutual = [moments[4], moments[5]];
        let duration = Dd::from_f64(self.shutter.effective_duration_s());
        let result = CoherentInterferometerExposureSample {
            arm_power_density_w: [rounded(arms[0])?, rounded(arms[1])?],
            output_power_density_w: [rounded(outputs[0])?, rounded(outputs[1])?],
            mutual_coherence_density_w: C64::new(rounded(mutual[0])?, rounded(mutual[1])?),
            fringe_visibility: visibility(arms, mutual)?,
            arm_energy_density_j: [energy(arms[0], duration)?, energy(arms[1], duration)?],
            output_energy_density_j: [energy(outputs[0], duration)?, energy(outputs[1], duration)?],
            pupil_evaluations: evaluations,
            pair_contractions: contractions,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    // Original DD moments for point and finite-area receivers. Relay offsets
    // stay separate until field evaluation; propagation retains their DD sum.
    fn interferometer_moments(
        &self,
        receiver: CoherentDiffractionInterferometer,
        relay_offset: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Dd; 6], CoherentDiffractionError> {
        let amplitudes = receiver.arm_amplitudes();
        let active = amplitudes.iter().filter(|&&a| a != Dd::ZERO).count();
        let mut arms = [Dd::ZERO; 2];
        let mut outputs = [Dd::ZERO; 2];
        let mut mutual = [Dd::ZERO; 2];
        if self.waves != 0 && active != 0 {
            let coordinates = receiver.coordinates().map(|center| {
                std::array::from_fn(|axis| {
                    (Dd::from_f64(center[axis]) + Dd::from_f64(relay_offset[axis])).to_f64()
                })
            });
            let rotations =
                self.interferometer_rotations_offset(receiver, relay_offset, cancelled)?;
            let identity = [[Dd::ONE, Dd::ZERO]; MAX_COHERENT_DIFFRACTION_TONES];
            let mixing_amplitude = Dd::from_f64(0.5).sqrt();
            for bank in 0..=self.factors {
                poll(cancelled)?;
                let mut fields = [[[Dd::ZERO; 2]; MAX_COHERENT_DIFFRACTION_TONES]; 2];
                for (tone, &group) in self.tones.iter().zip(&self.tone_groups) {
                    poll(cancelled)?;
                    let waves = if bank == 0 {
                        tone.illumination.mean
                    } else {
                        tone.illumination.factors[bank - 1].waves
                    };
                    for arm in 0..2 {
                        if amplitudes[arm] != Dd::ZERO {
                            let field = tone.field(waves, coordinates[arm], cancelled)?;
                            for axis in 0..2 {
                                fields[arm][group][axis] =
                                    fields[arm][group][axis] + field[axis] * amplitudes[arm];
                            }
                        }
                    }
                }
                let mut bank_arms = [Dd::ZERO; 2];
                for arm in 0..2 {
                    if amplitudes[arm] != Dd::ZERO {
                        for i in 0..self.groups.len() {
                            poll(cancelled)?;
                            fields[arm][i] = mul(fields[arm][i], rotations[arm][i]);
                        }
                        let mut component = fields[arm];
                        bank_arms[arm] =
                            self.contract_component(&mut component, &identity, cancelled)?;
                        arms[arm] = arms[arm] + bank_arms[arm];
                    }
                }
                if active == 1 {
                    let power = (bank_arms[0] + bank_arms[1]) * Dd::from_f64(0.5);
                    outputs[0] = outputs[0] + power;
                    outputs[1] = outputs[1] + power;
                    continue;
                }
                for (port, sign) in [Dd::ONE, -Dd::ONE].into_iter().enumerate() {
                    let mut component = [[Dd::ZERO; 2]; MAX_COHERENT_DIFFRACTION_TONES];
                    for i in 0..self.groups.len() {
                        for axis in 0..2 {
                            // Normalize the field BEFORE squaring, so a finite
                            // bright output cannot overflow at an unscaled port.
                            component[i][axis] =
                                (fields[0][i][axis] + sign * fields[1][i][axis]) * mixing_amplitude;
                        }
                    }
                    outputs[port] = outputs[port]
                        + self.contract_component(&mut component, &identity, cancelled)?;
                }
                let cross = self.interferometer_cross(&fields[0], &fields[1], cancelled)?;
                // Refuse a failed bank rather than hiding it under another
                // independent source's positive power.
                visibility(bank_arms, cross)?;
                for axis in 0..2 {
                    mutual[axis] = mutual[axis] + cross[axis];
                }
            }
        }
        Ok([
            arms[0], arms[1], outputs[0], outputs[1], mutual[0], mutual[1],
        ])
    }

    fn interferometer_work(
        &self,
        active: usize,
        maximum: usize,
    ) -> Result<usize, CoherentDiffractionError> {
        let quantity = "interferometer exposure pair contractions";
        let work = if active == 0 || self.waves == 0 {
            0
        } else {
            let per_bank = if active == 1 {
                self.power_component_work()
            } else {
                self.power_component_work()
                    .checked_mul(4)
                    .and_then(|x| x.checked_add(self.cross_component_work()))
                    .ok_or_else(|| limit(quantity, maximum))?
            };
            per_bank
                .checked_mul(
                    self.factors
                        .checked_add(1)
                        .ok_or_else(|| limit(quantity, maximum))?,
                )
                .ok_or_else(|| limit(quantity, maximum))?
        };
        admit_count(quantity, work, maximum)?;
        Ok(work)
    }

    pub(super) fn interferometer_rotations(
        &self,
        receiver: CoherentDiffractionInterferometer,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Fields; 2], CoherentDiffractionError> {
        self.interferometer_rotations_offset(receiver, [0.0; 2], cancelled)
    }

    pub(super) fn interferometer_rotations_offset(
        &self,
        receiver: CoherentDiffractionInterferometer,
        relay_offset: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Fields; 2], CoherentDiffractionError> {
        let centers = receiver.coordinates();
        let coordinates = centers.map(|center| {
            std::array::from_fn::<_, 2, _>(|axis| {
                Dd::from_f64(center[axis]) + Dd::from_f64(relay_offset[axis])
            })
        });
        // Drop only ONE common carrier/coordinate phase from both arms.
        // Independently dropping each arm's reference chirp would erase the
        // monochromatic spatial interference this receiver is measuring.
        let distance = Dd::from_f64(self.sensor_distance_m)
            * (Dd::ONE
                + (coordinates[0][0] * coordinates[0][0] + coordinates[0][1] * coordinates[0][1])
                    * Dd::from_f64(0.5));
        rounded(distance)?;
        // The shared offset cancels from the coordinate difference. Factoring
        // |c1+u|²-|c0+u|² avoids subtracting two independently rounded chirps.
        let radius_difference = (0..2).fold(Dd::ZERO, |sum, axis| {
            sum + (Dd::from_f64(centers[1][axis]) - Dd::from_f64(centers[0][axis]))
                * (coordinates[1][axis] + coordinates[0][axis])
        });
        let path = Dd::from_f64(self.sensor_distance_m) * radius_difference * Dd::from_f64(0.5);
        rounded(path)?;
        let (sine, cosine) = det::sin_cos_dd(Dd::from_f64(receiver.relative_phase_rad()));
        let interface = [cosine, sine];
        let mut first = [[Dd::ZERO; 2]; MAX_COHERENT_DIFFRACTION_TONES];
        let mut second = [[Dd::ZERO; 2]; MAX_COHERENT_DIFFRACTION_TONES];
        for (i, group) in self.groups.iter().enumerate() {
            poll(cancelled)?;
            let common = frequency_length_phase(group.relative_optical_frequency_hz, distance)?;
            let (sine, cosine) = det::sin_cos_dd(common);
            first[i] = [cosine, sine];
            let chirp = frequency_length_phase(group.optical_frequency_hz, path)?;
            let delay = frequency_length_phase(
                Dd::from_f64(group.frequency_hz),
                Dd::from_f64(receiver.optical_path_difference_m()),
            )?;
            let (cs, cc) = det::sin_cos_dd(chirp);
            let (ds, dc) = det::sin_cos_dd(delay);
            second[i] = mul(mul(mul(first[i], [cc, cs]), [dc, ds]), interface);
        }
        Ok([first, second])
    }

    pub(super) fn interferometer_cross(
        &self,
        left: &Fields,
        right: &Fields,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Dd; 2], CoherentDiffractionError> {
        poll(cancelled)?;
        if let Some(band) = &self.band {
            let count = self.groups.len();
            let cross = band.mean_cross_parts_dd(&left[..count], &right[..count])?;
            poll(cancelled)?;
            rounded(cross[0])?;
            rounded(cross[1])?;
            return Ok(cross);
        }
        let mut cross = [Dd::ZERO; 2];
        let two = Dd::from_f64(2.0);
        for i in 0..self.groups.len() {
            poll(cancelled)?;
            let diagonal = mul(left[i], [right[i][0], -right[i][1]]);
            cross[0] = cross[0] + diagonal[0];
            cross[1] = cross[1] + diagonal[1];
            for j in 0..i {
                poll(cancelled)?;
                let kernel = self.kernels[i * (i - 1) / 2 + j];
                // The native kernel returns half the real analytic product.
                // Rotating the right operand by +i extracts its imaginary
                // part; the reversed-frequency term is the conjugate pair.
                cross[0] = cross[0]
                    + two
                        * (kernel.mean_real_product_parts_dd(left[i], right[j])?
                            + kernel.mean_real_product_parts_dd(right[i], left[j])?);
                cross[1] = cross[1]
                    + two
                        * (kernel
                            .mean_real_product_parts_dd(left[i], [-right[j][1], right[j][0]])?
                            - kernel
                                .mean_real_product_parts_dd(right[i], [-left[j][1], left[j][0]])?);
            }
        }
        rounded(cross[0])?;
        rounded(cross[1])?;
        Ok(cross)
    }
}

// For exp(-i*2*pi*f*t), propagation through optical length L contributes
// exp(+i*2*pi*f*L/c). Reduce each exact scalar product modulo integer c before
// dividing or multiplying by TAU. A rounded c/f wavelength, a huge radian
// product, or collapsing the DD low components first can create false beats.
pub(super) fn frequency_length_phase(
    frequency_hz: Dd,
    optical_length_m: Dd,
) -> Result<Dd, CoherentDiffractionError> {
    let mut remainder = Dd::ZERO;
    for f in [frequency_hz.hi, frequency_hz.lo] {
        for length in [optical_length_m.hi, optical_length_m.lo] {
            let product = Dd::from_f64(f) * Dd::from_f64(length);
            rounded(product)?;
            remainder = remainder + Dd::from_pair(product.hi % C_M_S, product.lo % C_M_S);
        }
    }
    let phase = remainder / Dd::from_f64(C_M_S) * TAU;
    rounded(phase)?;
    Ok(phase)
}

pub(super) fn energy(power: Dd, duration: Dd) -> Result<f64, CoherentDiffractionError> {
    let value = power * duration;
    if power != Dd::ZERO && value == Dd::ZERO {
        return Err(CoherentDiffractionError::Unrepresentable);
    }
    rounded(value)
}
