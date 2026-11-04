//! Two-point electric coherence and ideal, explicitly authored field analyzers.
//!
//! Prepared Fourier coefficients already contain their inverse-FFT 1/N. Direct
//! synthesis uses the envelope Fourier phases and the authored transverse carrier,
//! with no further normalization. Coordinates observe the original periodic,
//! band-limited envelope and its carrier, not an isolated aperture. The existing
//! Maxwell mode owner supplies all longitudinal,
//! reflected, transmitted and material-dependent phases.

use super::*;

type ComplexDd = [Dd; 2];
type Electric = [ComplexDd; 3];
type Tensor = [[ComplexDd; 3]; 3];
const ZERO_COMPLEX: ComplexDd = [Dd::ZERO; 2];
const ZERO_TENSOR: Tensor = [[ZERO_COMPLEX; 3]; 3];

/// An exterior sample in the original beam's Cartesian coordinates, in metres.
/// Incident/reflected/total samples lie at z=-distance_from_face_m; transmitted
/// samples lie beyond the last face. Distance must be finite and nonnegative;
/// x/y may be any finite coordinates in the envelope's continuation and carrier.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumObservation {
    pub side: AngularSpectrumSide,
    pub distance_from_face_m: f64,
    pub position_m: [f64; 2],
}

/// Query limits, checked before any modal reconstruction. Exactly 2*M*N modes
/// are reconstructed for M source factors and N original Fourier modes. No FFT,
/// material solve, heap allocation or dense spatial covariance is needed.
///
/// Numerical workspace counts the returned tensors, two synthesized vector
/// fields and one mode's E/H buffers. Scalar compiler temporaries, borrowed
/// beam/factor storage and previously retained results are excluded. Both
/// queries use the same conservative bound, including the smaller direct sum.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumCorrelationLimits {
    pub maximum_mode_evaluations: usize,
    pub maximum_workspace_bytes: usize,
}

/// Cycle/ensemble electric second moments between two observations, in V²/m².
/// Cross entry (i,j) is W_ab[i,j] = sum_m E_a,m[i] E_b,m[j]* / 2.
/// The reverse observation order gives the conjugate transpose. Local blocks
/// are Hermitian; the cross block generally is not. No mean field is invented.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumMutualCoherence {
    at_a: Tensor,
    at_b: Tensor,
    cross: Tensor,
}

impl AngularSpectrumMutualCoherence {
    pub fn electric_at_a_v2_m2(&self) -> [[C64; 3]; 3] {
        rounded_tensor(self.at_a)
    }
    pub fn electric_at_b_v2_m2(&self) -> [[C64; 3]; 3] {
        rounded_tensor(self.at_b)
    }
    pub fn electric_cross_v2_m2(&self) -> [[C64; 3]; 3] {
        rounded_tensor(self.cross)
    }

    /// Analyze a = w_a† E_a and b = w_b† E_b using explicit, finite complex
    /// Cartesian weights. Weights are dimensionless and need not be normalized;
    /// their conjugation is part of the analyzer, including its complex gain.
    ///
    /// The returned mutual correlation is Gamma_ab = <a b*>. Combining
    /// a + exp(+i phi)*b gives I_a + I_b + 2*Re(exp(-i phi)*Gamma_ab).
    /// For deep nulls use the source owner's direct recombination query, which
    /// combines each factor's fields before forming a positive square. Tensor
    /// contraction cannot guarantee relative accuracy arbitrarily near a null.
    pub fn project(
        &self,
        analyzers: [[C64; 3]; 2],
    ) -> Result<AngularSpectrumInterferenceStatistics, AngularSpectrumError> {
        admit_analyzers(analyzers)?;
        let a = contract(&self.at_a, analyzers[0], analyzers[0])[0];
        let b = contract(&self.at_b, analyzers[1], analyzers[1])[0];
        let cross = contract(&self.cross, analyzers[0], analyzers[1]);
        let arm_mean_square_v2_m2 = [nonnegative(a)?, nonnegative(b)?];
        let mutual_v2_m2 = rounded(cross)?;
        let [ia, ib] = arm_mean_square_v2_m2;
        let visibility = if ia == 0.0 || ib == 0.0 {
            if mutual_v2_m2 != C64::ZERO {
                return Err(bad("zero analyzer intensity has nonzero mutual coherence"));
            }
            if ia == 0.0 && ib == 0.0 {
                None
            } else {
                Some(0.0)
            }
        } else {
            // Cauchy-Schwarz is stronger than visibility <= 1 for unequal
            // arms. Divide by each root to avoid forming an overflowing Ia*Ib.
            let root_a = det::sqrt(ia);
            let root_b = det::sqrt(ib);
            let coherence = det::hypot(
                (mutual_v2_m2.re / root_a) / root_b,
                (mutual_v2_m2.im / root_a) / root_b,
            );
            if !coherence.is_finite() || coherence > 1.0 + 32.0 * f64::EPSILON {
                return Err(bad("projected mutual coherence violates Cauchy-Schwarz"));
            }
            // Scale before adding intensities or taking a complex norm.
            let scale = ia.max(ib);
            let value = 2.0 * det::hypot(mutual_v2_m2.re / scale, mutual_v2_m2.im / scale)
                / (ia / scale + ib / scale);
            if !value.is_finite() || value > 1.0 + 32.0 * f64::EPSILON {
                return Err(bad(
                    "projected mutual coherence violates its intensity bound",
                ));
            }
            Some(value.min(1.0))
        };
        Ok(AngularSpectrumInterferenceStatistics {
            arm_mean_square_v2_m2,
            mutual_v2_m2,
            visibility,
        })
    }
}

/// Electric analyzer moments, not irradiance, detector watts or a detector
/// absorption/back-action model. Visibility is 2*|Gamma_ab|/(I_a+I_b) for an
/// ideal relative-phase scan. It is None when both represented intensities are
/// zero, and zero when only one arm is dark. Final roundoff above one is bounded
/// and clipped; a larger Cauchy-Schwarz violation refuses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumInterferenceStatistics {
    pub arm_mean_square_v2_m2: [f64; 2],
    pub mutual_v2_m2: C64,
    pub visibility: Option<f64>,
}

impl AngularSpectrumCoherence<'_, '_, '_> {
    /// Exact number of mode reconstructions admitted by either two-point query.
    pub fn correlation_mode_evaluations(&self) -> Result<usize, AngularSpectrumError> {
        self.factors
            .len()
            .checked_mul(self.anchor.modes.len())
            .and_then(|work| work.checked_mul(2))
            .ok_or_else(|| bad("two-point coherence work overflow"))
    }

    /// Fixed logical numerical workspace; see AngularSpectrumCorrelationLimits.
    pub fn correlation_workspace_bytes(&self) -> usize {
        core::mem::size_of::<AngularSpectrumMutualCoherence>()
            + 2 * core::mem::size_of::<Electric>()
            + 2 * core::mem::size_of::<[C64; 3]>()
    }

    /// Propagate the original coherent factors to two physical points, then
    /// form their local and mutual electric moments. Different sides/depths
    /// keep the same original film response and material snapshot. Runtime is
    /// O(factors * modes), without reconstructing or retaining FFT planes.
    pub fn mutual_coherence(
        &self,
        a: AngularSpectrumObservation,
        b: AngularSpectrumObservation,
        limits: AngularSpectrumCorrelationLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumMutualCoherence, AngularSpectrumError> {
        self.admit_correlation([a, b], limits, &mut cancelled)?;
        let mut result = AngularSpectrumMutualCoherence {
            at_a: ZERO_TENSOR,
            at_b: ZERO_TENSOR,
            cross: ZERO_TENSOR,
        };
        for factor in self.factors {
            let [ea, eb] = synthesize_pair(factor, [a, b], &mut cancelled)?;
            accumulate_tensor(&mut result.at_a, ea, ea);
            accumulate_tensor(&mut result.at_b, eb, eb);
            accumulate_tensor(&mut result.cross, ea, eb);
        }
        for tensor in [&result.at_a, &result.at_b, &result.cross] {
            for row in tensor {
                for &entry in row {
                    rounded(entry)?;
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Directly evaluate <|w_a† E_a + exp(+i phi) w_b† E_b|²>/2, summing
    /// mutually uncorrelated factors only AFTER each positive square. The 1/2
    /// is the single peak-phasor cycle average, not another ensemble weight.
    ///
    /// Analyzer weights and relative_phase_rad must be finite. This preserves
    /// weak destructive residuals without subtracting rounded arm intensities
    /// or tensor entries. It retains the authored binary64 geometry, modal
    /// amplitudes and gains; it does not recover their physical uncertainty or
    /// guarantee relative accuracy at every cancellation. Units are V²/m².
    pub fn recombined_electric_mean_square_v2_m2(
        &self,
        a: AngularSpectrumObservation,
        b: AngularSpectrumObservation,
        analyzers: [[C64; 3]; 2],
        relative_phase_rad: f64,
        limits: AngularSpectrumCorrelationLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        self.admit_correlation([a, b], limits, &mut cancelled)?;
        admit_analyzers(analyzers)?;
        if !relative_phase_rad.is_finite() {
            return Err(bad("interferometer phase must be finite"));
        }
        let (sine, cosine) = det::sin_cos_dd(Dd::from_f64(relative_phase_rad));
        let phase = [cosine, sine];
        let mut sum = Dd::ZERO;
        for factor in self.factors {
            let [ea, eb] = synthesize_pair(factor, [a, b], &mut cancelled)?;
            let combined = plus(
                analyze(ea, analyzers[0]),
                multiply(phase, analyze(eb, analyzers[1])),
            );
            sum = sum + multiply(half(combined), conjugate(combined))[0];
        }
        let value = nonnegative(sum)?;
        poll(&mut cancelled)?;
        Ok(value)
    }

    fn admit_correlation(
        &self,
        observations: [AngularSpectrumObservation; 2],
        limits: AngularSpectrumCorrelationLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AngularSpectrumError> {
        poll(cancelled)?;
        // Darkness does not make invalid geometry or work limits admissible.
        for observation in observations {
            if !observation.distance_from_face_m.is_finite()
                || observation.distance_from_face_m < 0.0
                || observation.position_m.iter().any(|v| !v.is_finite())
            {
                return Err(bad(
                    "coherence observations need finite coordinates and nonnegative depth",
                ));
            }
        }
        if self.correlation_mode_evaluations()? > limits.maximum_mode_evaluations
            || self.correlation_workspace_bytes() > limits.maximum_workspace_bytes
        {
            return Err(bad(
                "two-point coherence exceeds its work or workspace budget",
            ));
        }
        Ok(())
    }
}

fn synthesize_pair(
    factor: &PlanarAngularSpectrum<'_, '_>,
    observations: [AngularSpectrumObservation; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Electric; 2], AngularSpectrumError> {
    poll(cancelled)?;
    let mut sums = [[ZERO_COMPLEX; 3]; 2];
    for (index, mode) in factor.modes.iter().enumerate() {
        if index % 32 == 0 {
            poll(cancelled)?;
        }
        let kx = factor.envelope_wavenumber(0, index % factor.grid.width);
        let ky = factor.envelope_wavenumber(1, index / factor.grid.width);
        for (sum, observation) in sums.iter_mut().zip(observations) {
            let (electric, _) =
                factor.mode_field(mode, observation.side, observation.distance_from_face_m)?;
            // Keep the original envelope offsets: adding a much larger carrier
            // in binary64 can erase distinct bins before phase reconstruction.
            let phase = Dd::from_f64(kx) * Dd::from_f64(observation.position_m[0])
                + Dd::from_f64(ky) * Dd::from_f64(observation.position_m[1]);
            if !phase.hi.is_finite() || !phase.lo.is_finite() {
                return Err(bad("two-point spatial phase is unrepresentable"));
            }
            let (sine, cosine) = det::sin_cos_dd(phase);
            for (value, amplitude) in sum.iter_mut().zip(electric) {
                *value = plus(*value, multiply(lift(amplitude), [cosine, -sine]));
            }
        }
    }
    // The common carrier multiplies each completed envelope once. It cancels
    // in local moments, but its relative phase survives between observations.
    if factor.carrier_rad_m != [0.0; 2] {
        for (field, observation) in sums.iter_mut().zip(observations) {
            let phase = Dd::from_f64(factor.carrier_rad_m[0])
                * Dd::from_f64(observation.position_m[0])
                + Dd::from_f64(factor.carrier_rad_m[1]) * Dd::from_f64(observation.position_m[1]);
            if !phase.hi.is_finite() || !phase.lo.is_finite() {
                return Err(bad("two-point carrier phase is unrepresentable"));
            }
            let (sine, cosine) = det::sin_cos_dd(phase);
            for value in field {
                *value = multiply(*value, [cosine, -sine]);
            }
        }
    }
    for field in sums {
        for value in field {
            rounded(value)?;
        }
    }
    poll(cancelled)?;
    Ok(sums)
}

fn accumulate_tensor(tensor: &mut Tensor, a: Electric, b: Electric) {
    for (row, ai) in tensor.iter_mut().zip(a) {
        for (entry, bj) in row.iter_mut().zip(b) {
            *entry = plus(*entry, multiply(half(ai), conjugate(bj)));
        }
    }
}

fn contract(tensor: &Tensor, a: [C64; 3], b: [C64; 3]) -> ComplexDd {
    let mut sum = ZERO_COMPLEX;
    for (row, ai) in tensor.iter().zip(a) {
        for (&entry, bj) in row.iter().zip(b) {
            sum = plus(
                sum,
                multiply(multiply(conjugate(lift(ai)), entry), lift(bj)),
            );
        }
    }
    sum
}

fn analyze(field: Electric, weights: [C64; 3]) -> ComplexDd {
    let mut sum = ZERO_COMPLEX;
    for (value, weight) in field.into_iter().zip(weights) {
        sum = plus(sum, multiply(conjugate(lift(weight)), value));
    }
    sum
}

fn admit_analyzers(analyzers: [[C64; 3]; 2]) -> Result<(), AngularSpectrumError> {
    if analyzers.iter().flatten().any(|&z| !finite(z)) {
        return Err(bad("electric analyzer weights must be finite"));
    }
    Ok(())
}

fn lift(z: C64) -> ComplexDd {
    [Dd::from_f64(z.re), Dd::from_f64(z.im)]
}
fn conjugate(z: ComplexDd) -> ComplexDd {
    [z[0], -z[1]]
}
fn plus(a: ComplexDd, b: ComplexDd) -> ComplexDd {
    [a[0] + b[0], a[1] + b[1]]
}
fn multiply(a: ComplexDd, b: ComplexDd) -> ComplexDd {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
fn half(z: ComplexDd) -> ComplexDd {
    z.map(|value| value * Dd::from_f64(0.5))
}
fn rounded(z: ComplexDd) -> Result<C64, AngularSpectrumError> {
    if z.iter().any(|v| !v.hi.is_finite() || !v.lo.is_finite()) {
        return Err(bad("two-point electric moment is unrepresentable"));
    }
    Ok(C64::new(z[0].to_f64(), z[1].to_f64()))
}
fn nonnegative(value: Dd) -> Result<f64, AngularSpectrumError> {
    rounded([value, Dd::ZERO])?;
    if value.lt(Dd::ZERO) {
        return Err(bad("electric mean square must remain nonnegative"));
    }
    Ok(value.to_f64())
}
fn rounded_tensor(tensor: Tensor) -> [[C64; 3]; 3] {
    tensor.map(|row| row.map(|z| C64::new(z[0].to_f64(), z[1].to_f64())))
}
