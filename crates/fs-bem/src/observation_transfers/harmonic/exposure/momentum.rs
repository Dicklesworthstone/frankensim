//! Finite-time wave-momentum transport on the original acoustic source clock.
//!
//! For a homogeneous, inviscid, nondispersive fluid at rest, the instantaneous
//! scalar-potential wave stress is
//! S = (p²/(2 rho c²) - rho u.u/2) I + rho u u.
//! In a source-free irrotational region it obeys
//! d_t(p u/c²) + div(S) = 0; symmetry also gives angular wave-momentum balance.
//! This is the wave pseudomomentum current. A finite-window surface flux is
//! NOT total mechanical traction, a body's impulse, or the complete momentum
//! of the fluid. Second-order pressure, mean flow, field storage, enclosure
//! and matching volume/surface quadrature are not inferred.
//! See Stone, Phys. Rev. E 62, 1341 (2000), especially sections IV and VII:
//! <https://arxiv.org/abs/cond-mat/9909315>.
//!
//! Pressure and all velocity components come from the original prepared
//! exposure. Spatial DD contraction precedes spectral storage; source
//! coherence, real-coordinate sum frequencies and proper-coordinate beat
//! frequencies retain their original meaning. No BEM solve or time sampling
//! is introduced, and queries allocate no memory.

use super::*;
use fs_fft::harmonic_window::{
    AngularHarmonicSpectrumTerm, HarmonicProductTone, HarmonicSpectrumTermLimits,
    MAX_PRODUCT_SPECTRUM_TERMS, PreparedHarmonicProductSpectrum, ProductSpectrumError,
};

mod stress;
use stress::spatial_pair;

/// Admission for a six-component wave-momentum surface or volume observation.
/// Counts include repeated/zero quadrature entries and silent spectral terms.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicMomentumLimits {
    pub maximum_quadrature_points: usize,
    /// Quadrature entries times sum of squared retained frequency counts.
    /// Each pair has a fixed six-channel spatial contraction.
    pub maximum_quadrature_frequency_pairs: usize,
    /// TOTAL scalar terms over all six channels, before coalescing: six per
    /// ordered frequency pair, with six more for each real-coordinate pair.
    pub maximum_generated_terms: usize,
    /// Conservative peak logical payload: copied geometry/factor laws,
    /// temporary angular terms, and all retained scalar spectrum capacities.
    /// Borrowed source fields, stack values and allocator overhead are excluded.
    pub maximum_workspace_bytes: usize,
}

/// Signed outgoing wave-momentum rates through the authored area vectors.
/// These are instantaneous or window-averaged pseudomomentum fluxes, not an
/// inferred force or torque on a body.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicWaveMomentumFlux {
    /// Linear wave momentum per second [N].
    pub linear_momentum_rate_n: [f64; 3],
    /// Angular wave momentum per second [N m], about the supplied origin.
    pub angular_momentum_rate_nm: [f64; 3],
}

/// Signed linear and angular wave momentum, retaining compensated low parts.
/// Surface integration yields transported momentum; volume integration yields
/// stored momentum. Neither quantity alone establishes a body's reaction.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicWaveMomentumTransfer {
    /// Linear wave momentum [N s].
    pub linear_momentum_ns: [Dd; 3],
    /// Angular wave momentum [N m s], about the declared observation origin.
    pub angular_momentum_nms: [Dd; 3],
}

/// Prepared ensemble-mean wave-momentum flux across an authored surface.
///
/// Six original scalar product-spectrum owners retain x/y/z linear flux and
/// x/y/z moment flux. Preparation costs O(Q sum_b F_b² + T log T) with a fixed
/// six-component factor, and O(Q + T + factors) logical memory. Queries use
/// only the coalesced tones, independently of the receiver count. There is no
/// dense covariance, receiver-by-tone array, or repeated field contraction.
/// Spatial/source uncertainty remains unresolved.
pub struct PreparedHarmonicSurfaceWaveMomentum<'a> {
    reference: &'a PreparedHarmonicObservation,
    source_clock: &'a str,
    quadrature: Vec<HarmonicSurfaceQuadrature>,
    moment_origin_m: [f64; 3],
    factor_statistics: Vec<HarmonicExposureFactorStatistics>,
    spectra: SixSpectra,
}

impl<'a> PreparedHarmonicSurfaceWaveMomentum<'a> {
    pub fn new(
        exposure: &PreparedHarmonicExposure<'a>,
        quadrature: &[HarmonicSurfaceQuadrature],
        moment_origin_m: [f64; 3],
        limits: HarmonicMomentumLimits,
    ) -> Result<Self, TransferError> {
        Self::new_cancellable(exposure, quadrature, moment_origin_m, limits, || false)
    }

    /// Admit before allocation. Refusal/cancellation leaves the exposure
    /// reusable and publishes no partial observer. Empty quadrature has the
    /// same original physical scope and an exactly zero flux.
    pub fn new_cancellable(
        exposure: &PreparedHarmonicExposure<'a>,
        quadrature: &[HarmonicSurfaceQuadrature],
        moment_origin_m: [f64; 3],
        limits: HarmonicMomentumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        poll(&mut cancelled)?;
        admit(quadrature.len(), MAX_SURFACE_QUADRATURE_POINTS)?;
        admit(quadrature.len(), limits.maximum_quadrature_points)?;
        if moment_origin_m.iter().any(|x| !x.is_finite()) {
            return Err(bad("acoustic wave-momentum origin is invalid"));
        }
        for sample in quadrature {
            poll(&mut cancelled)?;
            if sample.receiver >= exposure.reference.receivers.len()
                || sample.oriented_area_m2.iter().any(|a| !a.is_finite())
                || exposure.reference.receivers[sample.receiver]
                    .position_m
                    .iter()
                    .any(|x| !x.is_finite())
            {
                return Err(bad("acoustic wave-momentum quadrature is invalid"));
            }
        }
        let geometry_bytes = quadrature
            .len()
            .checked_mul(core::mem::size_of::<HarmonicSurfaceQuadrature>())
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
    pub fn quadrature(&self) -> &[HarmonicSurfaceQuadrature] {
        &self.quadrature
    }
    pub fn moment_origin_m(&self) -> [f64; 3] {
        self.moment_origin_m
    }
    pub fn factor_statistics(&self) -> &[HarmonicExposureFactorStatistics] {
        &self.factor_statistics
    }
    /// Six signed real observables: linear x/y/z [N], then angular x/y/z [N m].
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

    pub fn evaluate(&self, time_s: f64) -> Result<HarmonicWaveMomentumFlux, TransferError> {
        self.evaluate_cancellable(time_s, || false)
    }

    pub fn evaluate_cancellable(
        &self,
        time_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicWaveMomentumFlux, TransferError> {
        poll(&mut cancelled)?;
        let mut values = [0.0; 6];
        for (value, spectrum) in values.iter_mut().zip(&self.spectra.components) {
            *value = spectrum
                .evaluate_cancellable(time_s, &mut cancelled)
                .map_err(spectrum_error)?;
        }
        poll(&mut cancelled)?;
        Ok(flux(values))
    }

    /// Normalized detector mean. A shaped detector mean is not an unweighted
    /// physical-interval impulse; use transported_momentum_dd for that integral.
    pub fn mean_flux(
        &self,
        window: HarmonicWindow,
    ) -> Result<HarmonicWaveMomentumFlux, TransferError> {
        self.mean_flux_cancellable(window, || false)
    }

    pub fn mean_flux_cancellable(
        &self,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicWaveMomentumFlux, TransferError> {
        poll(&mut cancelled)?;
        let mut values = [0.0; 6];
        for (value, spectrum) in values.iter_mut().zip(&self.spectra.components) {
            *value = spectrum
                .mean_cancellable(window, &mut cancelled)
                .map_err(spectrum_error)?;
        }
        poll(&mut cancelled)?;
        Ok(flux(values))
    }

    /// Signed, unweighted wave momentum crossing the surface during a physical
    /// interval. Finite nonnegative duration is required; zero is exact zero.
    pub fn transported_momentum_dd(
        &self,
        start_s: f64,
        duration_s: f64,
    ) -> Result<HarmonicWaveMomentumTransfer, TransferError> {
        self.transported_momentum_dd_cancellable(start_s, duration_s, || false)
    }

    pub fn transported_momentum_dd_cancellable(
        &self,
        start_s: f64,
        duration_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicWaveMomentumTransfer, TransferError> {
        self.transported_momentum_from_clock_dd_cancellable(
            Dd::from_f64(start_s),
            duration_s,
            cancelled,
        )
    }

    /// Retain a clock epoch plus elapsed time without rounding their sum to f64.
    /// No emission onset, window shape or material reaction is inferred.
    pub fn transported_momentum_from_clock_dd(
        &self,
        start_s: Dd,
        duration_s: f64,
    ) -> Result<HarmonicWaveMomentumTransfer, TransferError> {
        self.transported_momentum_from_clock_dd_cancellable(start_s, duration_s, || false)
    }

    pub fn transported_momentum_from_clock_dd_cancellable(
        &self,
        start_s: Dd,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicWaveMomentumTransfer, TransferError> {
        poll(&mut cancelled)?;
        let mut values = [Dd::ZERO; 6];
        for (value, spectrum) in values.iter_mut().zip(&self.spectra.components) {
            *value = spectrum
                .causal_integral_from_clock_dd_cancellable(start_s, duration_s, 0.0, &mut cancelled)
                .map_err(spectrum_error)?
                .integral_dd();
        }
        poll(&mut cancelled)?;
        Ok(HarmonicWaveMomentumTransfer {
            linear_momentum_ns: [values[0], values[1], values[2]],
            angular_momentum_nms: [values[3], values[4], values[5]],
        })
    }
}

fn flux(v: [f64; 6]) -> HarmonicWaveMomentumFlux {
    HarmonicWaveMomentumFlux {
        linear_momentum_rate_n: [v[0], v[1], v[2]],
        angular_momentum_rate_nm: [v[3], v[4], v[5]],
    }
}

/// Shared six-channel contraction result. Each channel contains a complex
/// difference tone and, for a real factor, a complex sum tone. The caller
/// includes the real-product 1/2 factor in its spatial bilinear contraction.
pub(super) type WaveMomentumPair = ([[Dd; 2]; 6], [[Dd; 2]; 6]);

pub(super) struct SixSpectra {
    pub(super) components: [PreparedHarmonicProductSpectrum; 6],
    pub(super) generated_terms: usize,
    pub(super) quadrature_frequency_pairs: usize,
    pub(super) workspace_bytes: usize,
}

/// Reuse the original exposure's complete factor banks and scalar spectral
/// owner for surface flux or volume storage. Geometry is validated by the
/// caller; geometry_bytes is its copied logical payload. All allocation and
/// source/frequency work is admitted before the spatial callback runs.
pub(super) fn compile_six<C: FnMut() -> bool>(
    exposure: &PreparedHarmonicExposure<'_>,
    quadrature_len: usize,
    geometry_bytes: usize,
    limits: HarmonicMomentumLimits,
    cancelled: &mut C,
    mut spatial_pair: impl FnMut(&Line, &Line, bool, &mut C) -> Result<WaveMomentumPair, TransferError>,
) -> Result<SixSpectra, TransferError> {
    poll(cancelled)?;
    if !exposure.wave_momentum {
        return Err(bad(
            "acoustic wave momentum requires nondispersive original fields",
        ));
    }
    admit(quadrature_len, MAX_SURFACE_QUADRATURE_POINTS)?;
    admit(quadrature_len, limits.maximum_quadrature_points)?;
    let mut pairs = 0usize;
    let mut per_component = 0usize;
    if quadrature_len != 0 {
        for bank in core::iter::once(&exposure.mean).chain(exposure.factors.iter()) {
            poll(cancelled)?;
            let count = bank
                .lines
                .len()
                .checked_mul(bank.lines.len())
                .ok_or_else(overflow)?;
            pairs = pairs.checked_add(count).ok_or_else(overflow)?;
            let multiplicity = if bank.statistics == HarmonicExposureFactorStatistics::Real {
                2
            } else {
                1
            };
            per_component = per_component
                .checked_add(count.checked_mul(multiplicity).ok_or_else(overflow)?)
                .ok_or_else(overflow)?;
        }
    }
    let quadrature_frequency_pairs = quadrature_len.checked_mul(pairs).ok_or_else(overflow)?;
    let generated_terms = per_component.checked_mul(6).ok_or_else(overflow)?;
    admit(quadrature_frequency_pairs, MAX_PRODUCTS)?;
    admit(
        quadrature_frequency_pairs,
        limits.maximum_quadrature_frequency_pairs,
    )?;
    admit(generated_terms, MAX_PRODUCT_SPECTRUM_TERMS)?;
    admit(generated_terms, limits.maximum_generated_terms)?;
    let tone_bytes = per_component
        .checked_mul(core::mem::size_of::<HarmonicProductTone>())
        .ok_or_else(overflow)?;
    let workspace_bytes = generated_terms
        .checked_mul(
            core::mem::size_of::<AngularHarmonicSpectrumTerm>()
                + core::mem::size_of::<HarmonicProductTone>(),
        )
        .and_then(|n| n.checked_add(geometry_bytes))
        .and_then(|n| {
            n.checked_add(
                exposure
                    .factors
                    .len()
                    .checked_mul(core::mem::size_of::<HarmonicExposureFactorStatistics>())?,
            )
        })
        .ok_or_else(overflow)?;
    admit(workspace_bytes, limits.maximum_workspace_bytes)?;
    let mut terms = [
        storage(per_component)?,
        storage(per_component)?,
        storage(per_component)?,
        storage(per_component)?,
        storage(per_component)?,
        storage(per_component)?,
    ];
    if quadrature_len != 0 {
        for bank in core::iter::once(&exposure.mean).chain(exposure.factors.iter()) {
            poll(cancelled)?;
            let real = bank.statistics == HarmonicExposureFactorStatistics::Real;
            for a in &bank.lines {
                for b in &bank.lines {
                    poll(cancelled)?;
                    let (difference, sum) = spatial_pair(a, b, real, cancelled)?;
                    let left = Dd::from_f64(a.omega);
                    let right = Dd::from_f64(b.omega);
                    for c in 0..6 {
                        terms[c].push(AngularHarmonicSpectrumTerm {
                            angular_frequency_rad_s: left - right,
                            peak_amplitude: difference[c],
                        });
                        if real {
                            terms[c].push(AngularHarmonicSpectrumTerm {
                                angular_frequency_rad_s: left + right,
                                peak_amplitude: sum[c],
                            });
                        }
                    }
                }
            }
        }
    }
    let mut build = |terms| {
        PreparedHarmonicProductSpectrum::from_angular_terms_cancellable(
            terms,
            HarmonicSpectrumTermLimits {
                maximum_terms: per_component,
                maximum_workspace_bytes: tone_bytes,
            },
            &mut *cancelled,
        )
        .map_err(spectrum_error)
    };
    let [x, y, z, lx, ly, lz] = terms;
    let components = [
        build(x)?,
        build(y)?,
        build(z)?,
        build(lx)?,
        build(ly)?,
        build(lz)?,
    ];
    poll(cancelled)?;
    Ok(SixSpectra {
        components,
        generated_terms,
        quadrature_frequency_pairs,
        workspace_bytes,
    })
}

pub(super) fn spectrum_error(error: ProductSpectrumError) -> TransferError {
    match error {
        ProductSpectrumError::Cancelled => TransferError::Cancelled,
        ProductSpectrumError::Budget => {
            bad("acoustic momentum spectrum exceeded its admitted budget")
        }
        ProductSpectrumError::Allocation => bad("acoustic momentum spectrum allocation failed"),
        ProductSpectrumError::Window(_) => bad("acoustic momentum temporal quantity is invalid"),
    }
}

#[cfg(test)]
mod tests;
