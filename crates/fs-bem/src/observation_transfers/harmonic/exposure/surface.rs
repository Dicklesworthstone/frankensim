//! Oriented acoustic power through explicitly authored physical quadrature.
//!
//! Pressure and velocity remain the original BEM fields. The expected signed
//! power is sum_q area_q dot E[p_q(t) u_q(t)], in watts. Neither a plane-wave
//! impedance nor a closed surface, outward direction or absorption is inferred.
//! Spatial contraction precedes spectral storage, retaining DC, beats and the
//! real-coordinate carrier-sum terms without a receiver-by-tone expansion.

use super::*;
use fs_fft::harmonic_window::{
    AngularHarmonicSpectrumTerm, HarmonicProductTone, HarmonicSpectrumTermLimits,
    MAX_PRODUCT_SPECTRUM_TERMS, PreparedHarmonicProductSpectrum, ProductSpectrumError,
};

/// Fixed ceiling on authored quadrature entries, including zero/repeated entries.
pub const MAX_SURFACE_QUADRATURE_POINTS: usize = 65_536;

/// One explicitly oriented area weight on an ORIGINAL physical receiver.
///
/// The vector includes both area [m²] and orientation. Every component must be
/// finite; signed and zero vectors remain authored. Repeated receiver indices
/// are allowed for subdivision or separate oriented patches. No normal, area
/// normalization, geometry interpolation or closed-surface inference occurs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicSurfaceQuadrature {
    pub receiver: usize,
    pub oriented_area_m2: [f64; 3],
}

/// Whole preparation limits. Queries are allocation-free in the prepared tones.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicSurfacePowerLimits {
    pub maximum_quadrature_points: usize,
    /// Quadrature count times sum of squared retained line counts over all
    /// banks. Each entry has a fixed three-component DD p/u contraction.
    pub maximum_quadrature_frequency_pairs: usize,
    /// One difference tone per ordered frequency pair; one additional sum tone
    /// for a deterministic/real bank. Counts precede coalescing, include silent
    /// terms, and are independent of the number of quadrature points.
    pub maximum_generated_terms: usize,
    /// Conservative peak logical payload: copied quadrature/factor laws,
    /// temporary angular terms and the retained spectrum's admitted capacity.
    /// Borrowed physical fields, stack values and allocator overhead are excluded.
    pub maximum_workspace_bytes: usize,
}

/// Prepared expected acoustic power waveform across an authored oriented surface.
///
/// Each source factor spans its original sources, receivers and frequencies.
/// Real/proper laws are preserved through the physical p/u contraction. This is
/// an ensemble-mean power waveform, not a random realization or a fourth-moment
/// power-noise spectrum. Negative flux is retained. Flux alone does not identify
/// dissipated heat, net radiated power, or a body's mechanical power balance.
///
/// Preparation uses O(Q * sum_b(F_b²) + T log T) work and O(Q + T + factors)
/// memory, for Q samples, F_b frequencies per bank and T generated terms.
/// Fields are contracted spatially in DD before emitting a frequency-pair term.
/// No additional BEM solve, covariance matrix, sampling or propagation model is
/// introduced. Original source/solver and quadrature uncertainty is unresolved.
pub struct PreparedHarmonicSurfacePower<'a> {
    reference: &'a PreparedHarmonicObservation,
    source_clock: &'a str,
    quadrature: Vec<HarmonicSurfaceQuadrature>,
    factor_statistics: Vec<HarmonicExposureFactorStatistics>,
    spectrum: PreparedHarmonicProductSpectrum,
    quadrature_frequency_pairs: usize,
    generated_terms: usize,
    workspace_bytes: usize,
}

impl<'a> PreparedHarmonicSurfacePower<'a> {
    pub fn new(
        exposure: &PreparedHarmonicExposure<'a>,
        quadrature: &[HarmonicSurfaceQuadrature],
        limits: HarmonicSurfacePowerLimits,
    ) -> Result<Self, TransferError> {
        Self::new_cancellable(exposure, quadrature, limits, || false)
    }

    /// Admission precedes allocation. Cancellation/refusal publishes no partial
    /// surface. Empty quadrature is an explicit zero flux with original scope.
    #[allow(clippy::too_many_lines)]
    pub fn new_cancellable(
        exposure: &PreparedHarmonicExposure<'a>,
        quadrature: &[HarmonicSurfaceQuadrature],
        limits: HarmonicSurfacePowerLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        poll(&mut cancelled)?;
        admit(quadrature.len(), MAX_SURFACE_QUADRATURE_POINTS)?;
        admit(quadrature.len(), limits.maximum_quadrature_points)?;
        for sample in quadrature {
            poll(&mut cancelled)?;
            if sample.receiver >= exposure.reference.receivers.len()
                || sample.oriented_area_m2.iter().any(|v| !v.is_finite())
            {
                return Err(bad(
                    "BEM surface quadrature receiver or oriented area is invalid",
                ));
            }
        }
        let mut frequency_pairs = 0usize;
        let mut generated_terms = 0usize;
        if !quadrature.is_empty() {
            for bank in core::iter::once(&exposure.mean).chain(exposure.factors.iter()) {
                poll(&mut cancelled)?;
                let pairs = bank
                    .lines
                    .len()
                    .checked_mul(bank.lines.len())
                    .ok_or_else(overflow)?;
                frequency_pairs = frequency_pairs.checked_add(pairs).ok_or_else(overflow)?;
                let multiplicity = match bank.statistics {
                    HarmonicExposureFactorStatistics::Real => 2,
                    HarmonicExposureFactorStatistics::ProperComplex => 1,
                };
                generated_terms = generated_terms
                    .checked_add(pairs.checked_mul(multiplicity).ok_or_else(overflow)?)
                    .ok_or_else(overflow)?;
            }
        }
        let quadrature_frequency_pairs = quadrature
            .len()
            .checked_mul(frequency_pairs)
            .ok_or_else(overflow)?;
        admit(quadrature_frequency_pairs, MAX_PRODUCTS)?;
        admit(
            quadrature_frequency_pairs,
            limits.maximum_quadrature_frequency_pairs,
        )?;
        admit(generated_terms, MAX_PRODUCT_SPECTRUM_TERMS)?;
        admit(generated_terms, limits.maximum_generated_terms)?;
        let tone_bytes = generated_terms
            .checked_mul(core::mem::size_of::<HarmonicProductTone>())
            .ok_or_else(overflow)?;
        let workspace_bytes = quadrature
            .len()
            .checked_mul(core::mem::size_of::<HarmonicSurfaceQuadrature>())
            .and_then(|bytes| {
                bytes.checked_add(
                    exposure
                        .factors
                        .len()
                        .checked_mul(core::mem::size_of::<HarmonicExposureFactorStatistics>())?,
                )
            })
            .and_then(|bytes| {
                bytes.checked_add(
                    generated_terms
                        .checked_mul(core::mem::size_of::<AngularHarmonicSpectrumTerm>())?,
                )
            })
            .and_then(|bytes| bytes.checked_add(tone_bytes))
            .ok_or_else(overflow)?;
        admit(workspace_bytes, limits.maximum_workspace_bytes)?;

        let mut owned_quadrature = storage(quadrature.len())?;
        owned_quadrature.extend_from_slice(quadrature);
        let mut factor_statistics = storage(exposure.factors.len())?;
        factor_statistics.extend(exposure.factor_statistics());
        let mut terms = storage(generated_terms)?;
        if !quadrature.is_empty() {
            for bank in core::iter::once(&exposure.mean).chain(exposure.factors.iter()) {
                poll(&mut cancelled)?;
                for a in &bank.lines {
                    for b in &bank.lines {
                        poll(&mut cancelled)?;
                        let real = bank.statistics == HarmonicExposureFactorStatistics::Real;
                        let (difference, sum) =
                            spatial_pair(a, b, quadrature, real, &mut cancelled)?;
                        let left = Dd::from_f64(a.omega);
                        let right = Dd::from_f64(b.omega);
                        terms.push(AngularHarmonicSpectrumTerm {
                            angular_frequency_rad_s: left - right,
                            peak_amplitude: difference,
                        });
                        if real {
                            terms.push(AngularHarmonicSpectrumTerm {
                                angular_frequency_rad_s: left + right,
                                peak_amplitude: sum,
                            });
                        }
                    }
                }
            }
        }
        let spectrum = PreparedHarmonicProductSpectrum::from_angular_terms_cancellable(
            terms,
            HarmonicSpectrumTermLimits {
                maximum_terms: generated_terms,
                maximum_workspace_bytes: tone_bytes,
            },
            &mut cancelled,
        )
        .map_err(surface_error)?;
        poll(&mut cancelled)?;
        Ok(Self {
            reference: exposure.reference,
            source_clock: exposure.source_clock,
            quadrature: owned_quadrature,
            factor_statistics,
            spectrum,
            quadrature_frequency_pairs,
            generated_terms,
            workspace_bytes,
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
        ProducerUncertainty::Unresolved {
            source_hash: self.reference.producer_hash,
        }
    }
    pub fn quadrature(&self) -> &[HarmonicSurfaceQuadrature] {
        &self.quadrature
    }
    pub fn factor_statistics(&self) -> &[HarmonicExposureFactorStatistics] {
        &self.factor_statistics
    }
    /// Real observable in watts on the original common source clock.
    pub fn spectrum(&self) -> &PreparedHarmonicProductSpectrum {
        &self.spectrum
    }
    pub fn quadrature_frequency_pairs(&self) -> usize {
        self.quadrature_frequency_pairs
    }
    pub fn generated_terms(&self) -> usize {
        self.generated_terms
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Expected signed instantaneous power [W], retaining all admitted beats.
    pub fn evaluate(&self, time_s: f64) -> Result<f64, TransferError> {
        self.evaluate_cancellable(time_s, || false)
    }

    pub fn evaluate_cancellable(
        &self,
        time_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, TransferError> {
        self.spectrum
            .evaluate_cancellable(time_s, cancelled)
            .map_err(surface_error)
    }

    /// Normalized detector average [W]. A shaped detector is not a physical
    /// interval-energy integral; use energy_j_dd for unweighted transported energy.
    pub fn mean_power_w(&self, window: HarmonicWindow) -> Result<f64, TransferError> {
        self.mean_power_w_cancellable(window, || false)
    }

    pub fn mean_power_w_cancellable(
        &self,
        window: HarmonicWindow,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, TransferError> {
        self.spectrum
            .mean_cancellable(window, cancelled)
            .map_err(surface_error)
    }

    /// Signed unweighted transported energy [J] over a physical interval.
    /// Duration is finite and nonnegative; zero duration returns exact zero.
    pub fn energy_j_dd(&self, start_s: f64, duration_s: f64) -> Result<Dd, TransferError> {
        self.energy_j_dd_cancellable(start_s, duration_s, || false)
    }

    pub fn energy_j_dd_cancellable(
        &self,
        start_s: f64,
        duration_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Dd, TransferError> {
        self.energy_j_dd_from_clock_cancellable(Dd::from_f64(start_s), duration_s, cancelled)
    }

    /// Preserve a source-clock epoch plus elapsed time without rounding their
    /// sum to f64. This does not infer an onset or alter the authored phase clock.
    pub fn energy_j_dd_from_clock(
        &self,
        start_s: Dd,
        duration_s: f64,
    ) -> Result<Dd, TransferError> {
        self.energy_j_dd_from_clock_cancellable(start_s, duration_s, || false)
    }

    pub fn energy_j_dd_from_clock_cancellable(
        &self,
        start_s: Dd,
        duration_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Dd, TransferError> {
        self.spectrum
            .causal_integral_from_clock_dd_cancellable(start_s, duration_s, 0.0, cancelled)
            .map(|result| result.integral_dd())
            .map_err(surface_error)
    }
}

fn spatial_pair(
    a: &Line,
    b: &Line,
    quadrature: &[HarmonicSurfaceQuadrature],
    real: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<([Dd; 2], [Dd; 2]), TransferError> {
    let mut difference = [Dd::ZERO; 2];
    let mut sum = [Dd::ZERO; 2];
    for sample in quadrature {
        poll(cancelled)?;
        let p = component(&a.fields[sample.receiver], 0);
        let velocity = &b.fields[sample.receiver];
        let mut projected = [Dd::ZERO; 2];
        for c in 0..3 {
            let area = Dd::from_f64(sample.oriented_area_m2[c]);
            let u = component(velocity, c + 1);
            projected[0] = projected[0] + area * u[0];
            projected[1] = projected[1] + area * u[1];
        }
        let rr = p[0] * projected[0];
        let ii = p[1] * projected[1];
        let ri = p[0] * projected[1];
        let ir = p[1] * projected[0];
        difference[0] = difference[0] + rr + ii;
        difference[1] = difference[1] + ir - ri;
        if real {
            sum[0] = sum[0] + rr - ii;
            sum[1] = sum[1] + ri + ir;
        }
    }
    let half = Dd::from_f64(0.5);
    for value in difference.iter_mut().chain(&mut sum) {
        *value = *value * half;
        if !finite(*value) {
            return Err(bad("BEM surface p/u contraction is unrepresentable"));
        }
    }
    Ok((difference, sum))
}

fn surface_error(error: ProductSpectrumError) -> TransferError {
    match error {
        ProductSpectrumError::Cancelled => TransferError::Cancelled,
        ProductSpectrumError::Budget => bad("BEM surface spectrum exceeded its admitted budget"),
        ProductSpectrumError::Allocation => bad("BEM surface spectrum allocation failed"),
        ProductSpectrumError::Window(_) => bad("BEM surface spectrum temporal quantity is invalid"),
    }
}

#[cfg(test)]
mod tests;
