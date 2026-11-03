//! Finite-volume storage of the original acoustic wave momentum.
//!
//! In the homogeneous inviscid nondispersive rest medium, wave momentum density
//! is p*u/c² and its orbital moment is (x-origin) cross (p*u/c²). Together with
//! the symmetric wave stress these obey source-free local conservation for
//! linear irrotational fields. This is wave pseudomomentum, not total fluid
//! momentum or a body's force. Neither a closed control volume, matching
//! surface quadrature, interior source terms nor second-order mean flow is
//! inferred. See Stone, Phys. Rev. E 62, 1341 (2000), sections IV-VII.

use super::momentum::{SixSpectra, WaveMomentumPair, compile_six, spectrum_error};
use super::*;
use fs_fft::harmonic_window::PreparedHarmonicProductSpectrum;

/// An authored nonnegative volume weight at an ORIGINAL field receiver.
/// Repeated receivers and zero weights are allowed and count against budgets.
/// Geometry coverage and cubature accuracy remain the caller's responsibility.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicVolumeQuadrature {
    pub receiver: usize,
    pub volume_m3: f64,
}

/// Instantaneous ensemble-mean wave momentum stored in the authored volume.
/// These are endpoint quantities, not normalized detector-window averages.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicWaveMomentumStorage {
    /// Linear wave momentum [N s].
    pub linear_momentum_ns: [Dd; 3],
    /// Angular wave momentum [N m s], about the observer's declared origin.
    pub angular_momentum_nms: [Dd; 3],
}

/// Prepared volume storage and exact analytic endpoint changes on one clock.
///
/// Spatial p/u contraction precedes six scalar spectrum compilations, using the
/// same real/proper factor laws and admitted source fields as surface flux.
/// Queries allocate no memory and do not revisit receiver/source products.
/// Cancellation or refusal never publishes a partial observer or six-vector.
pub struct PreparedHarmonicVolumeWaveMomentum<'a> {
    reference: &'a PreparedHarmonicObservation,
    source_clock: &'a str,
    quadrature: Vec<HarmonicVolumeQuadrature>,
    moment_origin_m: [f64; 3],
    factor_statistics: Vec<HarmonicExposureFactorStatistics>,
    spectra: SixSpectra,
}

impl<'a> PreparedHarmonicVolumeWaveMomentum<'a> {
    pub fn new(
        exposure: &PreparedHarmonicExposure<'a>,
        quadrature: &[HarmonicVolumeQuadrature],
        moment_origin_m: [f64; 3],
        limits: HarmonicMomentumLimits,
    ) -> Result<Self, TransferError> {
        Self::new_cancellable(exposure, quadrature, moment_origin_m, limits, || false)
    }

    pub fn new_cancellable(
        exposure: &PreparedHarmonicExposure<'a>,
        quadrature: &[HarmonicVolumeQuadrature],
        moment_origin_m: [f64; 3],
        limits: HarmonicMomentumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        poll(&mut cancelled)?;
        admit(quadrature.len(), MAX_SURFACE_QUADRATURE_POINTS)?;
        admit(quadrature.len(), limits.maximum_quadrature_points)?;
        if moment_origin_m.iter().any(|x| !x.is_finite()) {
            return Err(bad("acoustic volume momentum origin is invalid"));
        }
        let c = Dd::from_f64(exposure.reference.medium.sound_speed);
        let inverse_c2 = (Dd::ONE / c) / c;
        if !finite(inverse_c2) || inverse_c2 == Dd::ZERO {
            return Err(bad(
                "acoustic volume momentum medium scale is unrepresentable",
            ));
        }
        for sample in quadrature {
            poll(&mut cancelled)?;
            if sample.receiver >= exposure.reference.receivers.len()
                || !sample.volume_m3.is_finite()
                || sample.volume_m3 < 0.0
                || exposure.reference.receivers[sample.receiver]
                    .position_m
                    .iter()
                    .any(|x| !x.is_finite())
            {
                return Err(bad("acoustic volume momentum quadrature is invalid"));
            }
            let scale = Dd::from_f64(sample.volume_m3) * inverse_c2 * Dd::from_f64(0.5);
            if !finite(scale) || (sample.volume_m3 != 0.0 && scale == Dd::ZERO) {
                return Err(bad(
                    "acoustic volume momentum quadrature scale is unrepresentable",
                ));
            }
        }
        let geometry_bytes = quadrature
            .len()
            .checked_mul(core::mem::size_of::<HarmonicVolumeQuadrature>())
            .ok_or_else(overflow)?;
        let spectra = compile_six(
            exposure,
            quadrature.len(),
            geometry_bytes,
            limits,
            &mut cancelled,
            |a, b, real, cancelled| {
                spatial_pair(
                    a,
                    b,
                    exposure.reference,
                    quadrature,
                    moment_origin_m,
                    inverse_c2,
                    real,
                    cancelled,
                )
            },
        )?;
        let mut owned_quadrature = storage(quadrature.len())?;
        owned_quadrature.extend_from_slice(quadrature);
        let mut factor_statistics = storage(exposure.factors.len())?;
        factor_statistics.extend(exposure.factor_statistics());
        poll(&mut cancelled)?;
        Ok(Self {
            reference: exposure.reference,
            source_clock: exposure.source_clock,
            quadrature: owned_quadrature,
            moment_origin_m,
            factor_statistics,
            spectra,
        })
    }

    pub fn physical_receivers(&self) -> &[ReceiverSpec] {
        &self.reference.receivers
    }
    pub fn sources(&self) -> &[SourceId] {
        &self.reference.sources
    }
    pub fn source_clock(&self) -> &str {
        self.source_clock
    }
    pub fn producer_hash(&self) -> ContentHash {
        self.reference.producer_hash
    }
    pub fn producer_uncertainty(&self) -> ProducerUncertainty {
        self.reference.producer_uncertainty()
    }
    pub fn quadrature(&self) -> &[HarmonicVolumeQuadrature] {
        &self.quadrature
    }
    pub fn moment_origin_m(&self) -> [f64; 3] {
        self.moment_origin_m
    }
    pub fn factor_statistics(&self) -> &[HarmonicExposureFactorStatistics] {
        &self.factor_statistics
    }
    /// Linear x/y/z storage [N s], then angular x/y/z storage [N m s].
    pub fn spectra(&self) -> &[PreparedHarmonicProductSpectrum; 6] {
        &self.spectra.components
    }
    pub fn quadrature_frequency_pairs(&self) -> usize {
        self.spectra.quadrature_frequency_pairs
    }
    pub fn generated_terms(&self) -> usize {
        self.spectra.generated_terms
    }
    pub fn workspace_bytes(&self) -> usize {
        self.spectra.workspace_bytes
    }

    pub fn stored_momentum_dd(
        &self,
        time_s: f64,
    ) -> Result<HarmonicWaveMomentumStorage, TransferError> {
        self.stored_momentum_from_clock_dd(Dd::from_f64(time_s))
    }

    /// Preserve the source epoch plus elapsed time and the resulting low parts.
    pub fn stored_momentum_from_clock_dd(
        &self,
        time_s: Dd,
    ) -> Result<HarmonicWaveMomentumStorage, TransferError> {
        self.stored_momentum_from_clock_dd_cancellable(time_s, || false)
    }

    pub fn stored_momentum_from_clock_dd_cancellable(
        &self,
        time_s: Dd,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicWaveMomentumStorage, TransferError> {
        poll(&mut cancelled)?;
        let mut values = [Dd::ZERO; 6];
        for (value, spectrum) in values.iter_mut().zip(&self.spectra.components) {
            *value = spectrum
                .evaluate_from_clock_dd_cancellable(time_s, &mut cancelled)
                .map_err(spectrum_error)?;
        }
        poll(&mut cancelled)?;
        Ok(HarmonicWaveMomentumStorage {
            linear_momentum_ns: [values[0], values[1], values[2]],
            angular_momentum_nms: [values[3], values[4], values[5]],
        })
    }

    pub fn momentum_change_dd(
        &self,
        start_s: f64,
        duration_s: f64,
    ) -> Result<HarmonicWaveMomentumTransfer, TransferError> {
        self.momentum_change_from_clock_dd(Dd::from_f64(start_s), duration_s)
    }

    /// End minus start storage, calculated tone by tone without subtracting
    /// rounded large DC storage or collapsing a short interval into its epoch.
    /// Duration is finite and nonnegative; zero duration is exactly zero.
    /// Add this to outward transported wave momentum for an explicitly matched
    /// source-free control-volume balance. Matching quadrature and enclosure
    /// remain caller responsibilities; the result is not a body impulse.
    pub fn momentum_change_from_clock_dd(
        &self,
        start_s: Dd,
        duration_s: f64,
    ) -> Result<HarmonicWaveMomentumTransfer, TransferError> {
        self.momentum_change_from_clock_dd_cancellable(start_s, duration_s, || false)
    }

    pub fn momentum_change_from_clock_dd_cancellable(
        &self,
        start_s: Dd,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicWaveMomentumTransfer, TransferError> {
        poll(&mut cancelled)?;
        let mut values = [Dd::ZERO; 6];
        for (value, spectrum) in values.iter_mut().zip(&self.spectra.components) {
            *value = spectrum
                .increment_from_clock_dd_cancellable(start_s, duration_s, &mut cancelled)
                .map_err(spectrum_error)?;
        }
        poll(&mut cancelled)?;
        Ok(HarmonicWaveMomentumTransfer {
            linear_momentum_ns: [values[0], values[1], values[2]],
            angular_momentum_nms: [values[3], values[4], values[5]],
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn spatial_pair(
    a: &Line,
    b: &Line,
    reference: &PreparedHarmonicObservation,
    quadrature: &[HarmonicVolumeQuadrature],
    origin: [f64; 3],
    inverse_c2: Dd,
    real: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<WaveMomentumPair, TransferError> {
    let mut difference = [[Dd::ZERO; 2]; 6];
    let mut sum = [[Dd::ZERO; 2]; 6];
    for sample in quadrature {
        poll(cancelled)?;
        // Admission and work counts include zero entries. Do not form an
        // irrelevant huge lever arm or field product for an exact zero weight.
        if sample.volume_m3 == 0.0 {
            continue;
        }
        let p = component(&a.fields[sample.receiver], 0);
        let position = reference.receivers[sample.receiver].position_m;
        let arm: [Dd; 3] =
            core::array::from_fn(|i| Dd::from_f64(position[i]) - Dd::from_f64(origin[i]));
        let scale = Dd::from_f64(sample.volume_m3) * inverse_c2 * Dd::from_f64(0.5);
        let mut d = [[Dd::ZERO; 2]; 3];
        let mut s = [[Dd::ZERO; 2]; 3];
        for axis in 0..3 {
            let u = component(&b.fields[sample.receiver], axis + 1);
            d[axis] = [
                (p[0] * u[0] + p[1] * u[1]) * scale,
                (p[1] * u[0] - p[0] * u[1]) * scale,
            ];
            if real {
                s[axis] = [
                    (p[0] * u[0] - p[1] * u[1]) * scale,
                    (p[0] * u[1] + p[1] * u[0]) * scale,
                ];
            }
            for part in 0..2 {
                difference[axis][part] = difference[axis][part] + d[axis][part];
                sum[axis][part] = sum[axis][part] + s[axis][part];
            }
        }
        for axis in 0..3 {
            let j = (axis + 1) % 3;
            let k = (axis + 2) % 3;
            for part in 0..2 {
                difference[axis + 3][part] =
                    difference[axis + 3][part] + arm[j] * d[k][part] - arm[k] * d[j][part];
                if real {
                    sum[axis + 3][part] =
                        sum[axis + 3][part] + arm[j] * s[k][part] - arm[k] * s[j][part];
                }
            }
        }
    }
    if difference.iter().chain(&sum).flatten().any(|x| !finite(*x)) {
        return Err(bad(
            "acoustic volume p/u momentum contraction is unrepresentable",
        ));
    }
    Ok((difference, sum))
}

#[cfg(test)]
mod tests;
