//! Full spectral coherency through the original analytic shutter kernel.
//!
//! For each independent source bank, C_ab = <E_a(t) conj(E_b(t))>.
//! This is the average outer product, not the outer product of an averaged
//! field. The source clock, chromatic sensor phase and pupil transfer remain
//! owned by the existing Jones exposure. A finite shutter can depolarize
//! deterministic multicolor illumination by resolving its polarization beats.
//!
//! The two diagonal norms use the original Hermitian contraction. The complex
//! cross moment is contracted directly, never recovered by subtracting powers
//! measured through different analyzers. The cached kernel returns half of
//! Re[a conj(b) K_ij]; using i*b gives half of Im[a conj(b) K_ij]. For i>j the
//! cross entry is X_i conj(Y_j) K_ij + conj(Y_i conj(X_j) K_ij).
//! Thus both real terms add while their imaginary terms subtract.
//!
//! Camera right/up and incoming -forward define the same transverse frame as
//! monochromatic coherency, including native V = +2 Im(C_xy). The scalar
//! paraxial, analytic-field detector model and its no-claim boundaries remain.

use super::super::super::polarization::{
    PolarizedCoherentDiffractionPixelPower, PolarizedDiffractionCoherency, coherency, outer,
};
use super::{
    CoherentDiffractionError, CoherentDiffractionExposurePixelLimits, CoherentDiffractionPixel, Dd,
    JonesFields, MAX_COHERENT_DIFFRACTION_TONES, PreparedCoherentJonesDiffractionExposure,
    admit_coordinate, admit_count, limit, pixel, poll, rounded,
};

type Moments = [Dd; 4];

/// Shutter-averaged coherency per d(tan(theta_x)) d(tan(theta_y)).
///
/// Each result retains both axis powers, complex cross moment and native
/// Stokes. No single Jones field is assigned to a time-integrating detector.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedCoherentDiffractionExposureSample {
    /// Average outer product of the deterministic source mean, in W per area.
    pub mean_coherency_density_w: PolarizedDiffractionCoherency,
    /// Average contribution from independent covariance columns, in W per area.
    pub fluctuation_coherency_density_w: PolarizedDiffractionCoherency,
    /// Total ensemble/time-averaged coherency and Stokes density.
    pub coherency_density_w: PolarizedDiffractionCoherency,
}

/// Full polarization averaged over a physical pixel and uniform shutter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedCoherentDiffractionPixelExposure {
    /// Spatially integrated average coherency in watts, with exact work counts.
    pub average_power: PolarizedCoherentDiffractionPixelPower,
    /// Time-integrated deterministic-mean coherency in joules.
    pub mean_coherency_j: PolarizedDiffractionCoherency,
    /// Time-integrated covariance coherency in joules.
    pub fluctuation_coherency_j: PolarizedDiffractionCoherency,
    /// Total time-integrated coherency and signed Stokes energy in joules.
    pub coherency_j: PolarizedDiffractionCoherency,
    /// Real frequency-pair contractions over all source banks and pixel samples.
    pub pair_contractions: usize,
}

impl PreparedCoherentJonesDiffractionExposure<'_, '_> {
    /// Observe the full finite-shutter polarization without an analyzer.
    ///
    /// Equal frequencies coalesce before all products; independent source
    /// columns never acquire cross terms. Both complex off-diagonal terms of
    /// each frequency pair use the original analytic window and sensor phase.
    /// A common optical carrier is not time-stepped or rounded away.
    ///
    /// Full coherency costs n*(3*n+1) real contractions per source bank for n
    /// distinct frequencies: four diagonal moments and six per unequal pair.
    /// This is admitted against the original preparation's pair limit before
    /// pupil work. Detector-power queries retain their original smaller cost.
    pub fn observe_coherency(
        &self,
        at: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedCoherentDiffractionExposureSample, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        admit_coordinate(at)?;
        self.coherency_work()?;
        let [mean, fluctuations] = self.coherency_moments(at, &mut cancelled)?;
        let result = PolarizedCoherentDiffractionExposureSample {
            mean_coherency_density_w: physical(mean)?,
            fluctuation_coherency_density_w: physical(fluctuations)?,
            coherency_density_w: physical(std::array::from_fn(|i| mean[i] + fluctuations[i]))?,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Integrate full coherency over the original pixel and finite shutter.
    ///
    /// The existing positive spatial quadrature resolves both pupil bandwidth
    /// and chromatic sensor chirp. It integrates all matrix entries together,
    /// preserving temporal and spatial depolarization. Each vector wave evaluates
    /// its original pupil once per retained Jones-pupil mode and spatial sample.
    /// Unmasked waves have one mode. Pair work is
    /// preflighted against preparation and total pixel limits; queries allocate
    /// nothing and cancellation returns no partial result. Refine the existing
    /// spatial phase-span bound to assess quadrature error.
    pub fn integrate_pixel_coherency(
        &self,
        pixel: CoherentDiffractionPixel,
        limits: CoherentDiffractionExposurePixelLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedCoherentDiffractionPixelExposure, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        let contractions = self.coherency_work()?;
        let (subdivisions, quadrature_samples, pupil_evaluations, pair_contractions) = self
            .optics
            .pixel_work(pixel, limits, self.pupil_evaluations, contractions)?;
        let mut nonzero = [false; 4];
        let power: [Dd; 8] =
            pixel::integrate_moments(pixel, subdivisions, &mut cancelled, |at, cancelled| {
                let [mean, fluctuations] = self.coherency_moments(at, cancelled)?;
                for (flag, value) in
                    nonzero
                        .iter_mut()
                        .zip([mean[0], mean[1], fluctuations[0], fluctuations[1]])
                {
                    *flag |= !is_zero(value);
                }
                Ok(std::array::from_fn(|i| {
                    if i < 4 { mean[i] } else { fluctuations[i - 4] }
                }))
            })?;
        for (active, i) in nonzero.into_iter().zip([0, 1, 4, 5]) {
            if active && is_zero(power[i]) {
                return Err(CoherentDiffractionError::Unrepresentable);
            }
        }
        let mean = std::array::from_fn(|i| power[i]);
        let fluctuations = std::array::from_fn(|i| power[i + 4]);
        let total = std::array::from_fn(|i| mean[i] + fluctuations[i]);
        let duration = Dd::from_f64(self.shutter().duration_s());
        let result = PolarizedCoherentDiffractionPixelExposure {
            average_power: PolarizedCoherentDiffractionPixelPower {
                mean_coherency_w: physical(mean)?,
                fluctuation_coherency_w: physical(fluctuations)?,
                coherency_w: physical(total)?,
                subdivisions,
                quadrature_samples,
                pupil_evaluations,
            },
            mean_coherency_j: energy(mean, duration)?,
            fluctuation_coherency_j: energy(fluctuations, duration)?,
            coherency_j: energy(total, duration)?,
            pair_contractions,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn coherency_work(&self) -> Result<usize, CoherentDiffractionError> {
        let n = self.frequency_count();
        let contractions = (n * (3 * n + 1))
            .checked_mul(self.factors + 1)
            .ok_or_else(|| {
                limit(
                    "Jones coherency pair contractions",
                    self.maximum_pair_contractions,
                )
            })?;
        admit_count(
            "Jones coherency pair contractions",
            contractions,
            self.maximum_pair_contractions,
        )?;
        Ok(contractions)
    }

    fn coherency_moments(
        &self,
        at: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Moments; 2], CoherentDiffractionError> {
        let mut result = [[Dd::ZERO; 4]; 2];
        if self.pupil_evaluations == 0 {
            return Ok(result);
        }
        let rotations = self.optics.rotations(at, cancelled)?;
        for bank in 0..=self.factors {
            poll(cancelled)?;
            let mut fields = self.group_fields(bank, at, cancelled)?;
            let moment = self.contract_coherency(&mut fields, &rotations, cancelled)?;
            let index = usize::from(bank != 0);
            for i in 0..4 {
                result[index][i] = result[index][i] + moment[i];
            }
        }
        Ok(result)
    }

    pub(super) fn contract_coherency(
        &self,
        fields: &mut JonesFields,
        rotations: &[[Dd; 2]; MAX_COHERENT_DIFFRACTION_TONES],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Moments, CoherentDiffractionError> {
        let two = Dd::from_f64(2.0);
        let four = Dd::from_f64(4.0);
        let mut moment = [Dd::ZERO; 4];
        for i in 0..self.frequency_count() {
            poll(cancelled)?;
            let [cr, ci] = rotations[i];
            for field in &mut fields[i] {
                let [ar, ai] = *field;
                rounded(ar)?;
                rounded(ai)?;
                *field = [ar * cr - ai * ci, ar * ci + ai * cr];
                if field.iter().all(|&v| is_zero(v)) && (!is_zero(ar) || !is_zero(ai)) {
                    return Err(CoherentDiffractionError::Unrepresentable);
                }
            }
            let diagonal = outer(fields[i])?;
            for k in 0..4 {
                moment[k] = moment[k] + diagonal[k];
            }
            let [x, y] = fields[i];
            for j in 0..i {
                poll(cancelled)?;
                let kernel = self.optics.kernels[i * (i - 1) / 2 + j];
                let [u, v] = fields[j];
                moment[0] = moment[0] + four * kernel.mean_real_product_parts_dd(x, u)?;
                moment[1] = moment[1] + four * kernel.mean_real_product_parts_dd(y, v)?;
                moment[2] = moment[2]
                    + two
                        * (kernel.mean_real_product_parts_dd(x, v)?
                            + kernel.mean_real_product_parts_dd(y, u)?);
                moment[3] = moment[3]
                    + two
                        * (kernel.mean_real_product_parts_dd(x, [-v[1], v[0]])?
                            - kernel.mean_real_product_parts_dd(y, [-u[1], u[0]])?);
            }
        }
        // A failed bank cannot be hidden by adding another independent source.
        physical(moment)?;
        Ok(moment)
    }
}

fn is_zero(value: Dd) -> bool {
    value.hi == 0.0 && value.lo == 0.0
}

pub(super) fn energy(
    moment: Moments,
    duration: Dd,
) -> Result<PolarizedDiffractionCoherency, CoherentDiffractionError> {
    let energy = moment.map(|v| v * duration);
    for i in 0..2 {
        if !is_zero(moment[i]) && is_zero(energy[i]) {
            return Err(CoherentDiffractionError::Unrepresentable);
        }
    }
    physical(energy)
}

pub(super) fn physical(
    moment: Moments,
) -> Result<PolarizedDiffractionCoherency, CoherentDiffractionError> {
    let result = coherency(moment)?;
    let [x, y] = result.axis_powers();
    let cross = result.cross();
    // Rounded I/Q can hide a much weaker diagonal. Validate its correlation
    // directly as well, without squaring or multiplying extreme powers.
    if x == 0.0 || y == 0.0 {
        if cross.re != 0.0 || cross.im != 0.0 {
            return Err(CoherentDiffractionError::Unrepresentable);
        }
    } else {
        let large = x.max(y).sqrt();
        let small = x.min(y).sqrt();
        let correlation = fs_math::det::hypot(cross.re / large / small, cross.im / large / small);
        if correlation > 1.0 + 16.0 * f64::EPSILON {
            return Err(CoherentDiffractionError::Unrepresentable);
        }
    }
    Ok(result)
}

