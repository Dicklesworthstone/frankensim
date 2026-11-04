//! Measured broadband momentum on the original inviscid radiation-stress owner.

use super::super::welch::{peak, source_bin};
use super::*;
use fs_fft::welch::{WelchCrossSpectrum, WelchPeakBin};

/// One measured frequency bin and its original physical BEM transfer bank.
/// Ordered band entries must have strictly increasing bin indices. Missing
/// bins remain outside the requested band; they are never interpolated.
#[derive(Clone, Copy)]
pub struct HarmonicMomentumWelchLine<'a> {
    pub bin: usize,
    pub field: &'a PreparedHarmonicObservation,
}

/// Measured, band-integrated cycle-mean momentum flux [N, N m].
///
/// Every frequency uses its original p/u transfers and the same source kinds,
/// receiver geometry, fluid, boundary identity, surface and moment origin.
/// Measured source factors are propagated before the original radiation-stress
/// products; this retains spatial, incident/scattered and cross-source
/// interference. Frequencies are added as spectral powers. Welch does not
/// establish cross-frequency covariance, beat transients, finite-time impulse,
/// streaming, or statistical confidence.
///
/// Geometry and spectra are borrowed. Preparation and evaluation allocate no
/// heap memory. A fixed 64 KiB source buffer converts each segment coefficient
/// once, independently of the selected surface size. Only selected receivers
/// are evaluated; no dense source covariance or new BEM solve is formed.
/// Whole-band work is admitted once, including silent and repeated patches.
///
/// The original surface observer's physical scope applies. A closed outward
/// surface enclosing the complete system in a homogeneous nondispersive fluid
/// can establish its cycle-mean reaction as the negative outgoing flux. An
/// arbitrary patch does not establish a body force or a converged enclosure.
pub struct PreparedWelchSurfaceMomentum<'a> {
    spectrum: &'a WelchCrossSpectrum,
    lines: &'a [HarmonicMomentumWelchLine<'a>],
    quadrature: &'a [HarmonicSurfaceQuadrature],
    moment_origin_m: [f64; 3],
    products: usize,
}

impl<'a> PreparedWelchSurfaceMomentum<'a> {
    /// Bind an explicitly selected nonempty band to the original physical
    /// surface. DC, Nyquist, repeated bins, frequency mismatch, source/order
    /// substitution and dispersive/relaxing exterior stress are refused even
    /// for silent records. A record channel uses its source's original SI unit.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        spectrum: &'a WelchCrossSpectrum,
        lines: &'a [HarmonicMomentumWelchLine<'a>],
        recorded_order: &[SourceId],
        quadrature: &'a [HarmonicSurfaceQuadrature],
        moment_origin_m: [f64; 3],
        budget: HarmonicWorkBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        poll(&mut cancelled)?;
        cap(lines.len(), MAX_MODES)?;
        admit_products(0, budget)?;
        let reference = lines
            .first()
            .ok_or(TransferError::Invalid("measured momentum band is empty"))?
            .field;
        reference.momentum_surface(quadrature, moment_origin_m, &mut cancelled)?;
        let mut previous = None;
        let mut products = 0usize;
        for line in lines {
            poll(&mut cancelled)?;
            if previous.is_some_and(|bin| line.bin <= bin) {
                return Err(TransferError::Invalid(
                    "measured momentum bins must be strictly increasing",
                ));
            }
            previous = Some(line.bin);
            let field = line.field;
            if !field.radiation_momentum
                || field.producer_hash != reference.producer_hash
                || field.sources != reference.sources
                || field.receivers != reference.receivers
                || field.medium.density.to_bits() != reference.medium.density.to_bits()
                || field.medium.sound_speed.to_bits() != reference.medium.sound_speed.to_bits()
            {
                return Err(TransferError::Invalid(
                    "measured momentum band changed its original physical owner",
                ));
            }
            let view = source_bin(
                spectrum,
                line.bin,
                field.omega,
                &field.sources,
                recorded_order,
            )?;
            products = products
                .checked_add(work_products(
                    quadrature.len(),
                    field.sources.len(),
                    view.segment_count(),
                )?)
                .ok_or(TransferError::Invalid("measured momentum work overflow"))?;
            admit_products(products, budget)?;
        }
        poll(&mut cancelled)?;
        Ok(Self {
            spectrum,
            lines,
            quadrature,
            moment_origin_m,
            products,
        })
    }

    pub fn spectrum(&self) -> &WelchCrossSpectrum {
        self.spectrum
    }
    pub fn lines(&self) -> &[HarmonicMomentumWelchLine<'a>] {
        self.lines
    }
    pub fn quadrature(&self) -> &[HarmonicSurfaceQuadrature] {
        self.quadrature
    }
    pub fn moment_origin_m(&self) -> [f64; 3] {
        self.moment_origin_m
    }
    pub fn producer_hash(&self) -> ContentHash {
        self.lines[0].field.producer_hash
    }
    pub fn producer_uncertainty(&self) -> ProducerUncertainty {
        self.lines[0].field.producer_uncertainty()
    }
    pub fn work_products(&self) -> usize {
        self.products
    }

    pub fn evaluate(&self) -> Result<HarmonicMomentumFlux, TransferError> {
        self.evaluate_cancellable(|| false)
    }

    /// Accumulate the entire admitted band in compensated stress/torque sums,
    /// rounding only the final result. Output is already integrated over bin
    /// widths: do not multiply by df again. Cancellation publishes no partial
    /// force/torque and leaves this immutable observer reusable.
    pub fn evaluate_cancellable(
        &self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicMomentumFlux, TransferError> {
        poll(&mut cancelled)?;
        let mut sum = MomentumSum::default();
        let mut coefficients = [C64::ZERO; MAX_OBSERVATION_INPUTS];
        for line in self.lines {
            poll(&mut cancelled)?;
            let surface = PreparedHarmonicSurfaceMomentum {
                field: line.field,
                quadrature: self.quadrature,
                moment_origin_m: self.moment_origin_m,
            };
            let view = source_bin(
                self.spectrum,
                line.bin,
                line.field.omega,
                &line.field.sources,
                &line.field.sources,
            )?;
            accumulate(&surface, view, &mut coefficients, &mut sum, &mut cancelled)?;
        }
        poll(&mut cancelled)?;
        sum.finish()
    }
}

impl PreparedHarmonicSurfaceMomentum<'_> {
    /// Measured contribution of one positive interior frequency bin to
    /// cycle-mean radiation momentum. This is N/N m in the bin, not per Hz.
    /// The shared peak-bin conversion supplies sqrt(2 df) exactly once and
    /// the original stress supplies the peak-to-cycle-average factor of 1/2.
    /// Source phase and full near-field velocity correlations are retained.
    pub fn evaluate_welch_bin(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        recorded_order: &[SourceId],
        budget: HarmonicWorkBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicMomentumFlux, TransferError> {
        poll(&mut cancelled)?;
        let view = source_bin(
            spectrum,
            bin,
            self.field.omega,
            &self.field.sources,
            recorded_order,
        )?;
        admit_products(
            work_products(
                self.quadrature.len(),
                self.field.sources.len(),
                view.segment_count(),
            )?,
            budget,
        )?;
        let mut sum = MomentumSum::default();
        let mut coefficients = [C64::ZERO; MAX_OBSERVATION_INPUTS];
        accumulate(self, view, &mut coefficients, &mut sum, &mut cancelled)?;
        poll(&mut cancelled)?;
        sum.finish()
    }
}

fn accumulate(
    surface: &PreparedHarmonicSurfaceMomentum<'_>,
    view: WelchPeakBin<'_>,
    coefficients: &mut [C64; MAX_OBSERVATION_INPUTS],
    sum: &mut MomentumSum,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), TransferError> {
    let coefficients = &mut coefficients[..surface.field.sources.len()];
    for factor in view.factors() {
        poll(cancelled)?;
        for (i, (out, &coefficient)) in coefficients.iter_mut().zip(factor).enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            *out = peak(view, coefficient)?;
        }
        for sample in surface.quadrature {
            poll(cancelled)?;
            // The original compensated field reduction has a fixed 4096-source
            // ceiling. Pressure and velocity are summed before stress products.
            let field = surface.field.field(coefficients, sample.receiver)?;
            surface.add(sum, sample, field, 1.0)?;
        }
    }
    Ok(())
}

fn work_products(points: usize, sources: usize, segments: usize) -> Result<usize, TransferError> {
    // Two component scalings per source, then four transfer products per
    // source/receiver and seven independent pressure/velocity stress products.
    // Fixed 3x3 area/torque contractions accompany each admitted sample.
    sources
        .checked_mul(4)
        .and_then(|v| v.checked_add(7))
        .and_then(|v| v.checked_mul(points))
        .and_then(|v| v.checked_add(sources.checked_mul(2)?))
        .and_then(|v| v.checked_mul(segments))
        .ok_or(TransferError::Invalid("measured momentum work overflow"))
}

fn admit_products(products: usize, budget: HarmonicWorkBudget) -> Result<(), TransferError> {
    cap(budget.maximum_products, MAX_PRODUCTS)?;
    if budget.maximum_products == 0 {
        return Err(TransferError::Invalid("empty harmonic work allowance"));
    }
    cap(products, budget.maximum_products)
}

#[cfg(test)]
mod tests;
