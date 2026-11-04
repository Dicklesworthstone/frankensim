//! Measured source coherence through the original retained BEM transfers.
//!
//! A channel is the actual SI coordinate of its declared SourceId: structural
//! coordinate rate, aperture/point volume flow, or force on the fluid. The
//! caller supplies the recorded order explicitly; same-sized source banks are
//! not interchangeable. All channels share one physical clock and frame.

use super::*;
use fs_fft::welch::{WelchCrossSpectrum, WelchPeakBin};

pub(in crate::observation_transfers) fn source_bin<'a>(
    spectrum: &'a WelchCrossSpectrum,
    bin: usize,
    omega: f64,
    sources: &[SourceId],
    recorded_order: &[SourceId],
) -> Result<WelchPeakBin<'a>, TransferError> {
    let view = spectrum
        .peak_bin_rad_s(bin, omega)
        .map_err(|_| TransferError::Invalid("recorded source bin/frequency is not admissible"))?;
    if sources != recorded_order || view.channel_count() != sources.len() {
        return Err(TransferError::Invalid(
            "recorded source kinds and order must match the retained physical operator",
        ));
    }
    Ok(view)
}

pub(in crate::observation_transfers) fn peak(
    view: WelchPeakBin<'_>,
    coefficient: fs_fft::C64,
) -> Result<C64, TransferError> {
    let value = view
        .peak_factor(coefficient)
        .map_err(|_| TransferError::Invalid("recorded source peak factor is unrepresentable"))?;
    Ok(C64::new(value.re, value.im))
}

impl PreparedHarmonicObservation {
    /// Observe one measured source-spectrum bin through the original p/u bank.
    ///
    /// `recorded_order` identifies each history channel in its original SI
    /// coordinate; it must exactly equal `sources()`, including source kinds.
    /// The original source normalization, fluid, boundary geometry, scattering
    /// and molecular response remain those of this prepared owner. No new
    /// boundary solve, receiver quadrature or dense source covariance is formed.
    ///
    /// Results are this bin's pressure/velocity powers and signed intensity,
    /// not densities per hertz. The shared peak-bin view supplies sqrt(2*df)
    /// once; the existing statistics then apply the peak factor of one half.
    /// DC/Nyquist refuse. Segment factors need not be independent realizations;
    /// finite-record resolution and leakage remain, with no statistical
    /// confidence or cross-frequency covariance inferred.
    ///
    /// Allocation-free O(segments * sources * receivers), with explicit work
    /// admission and bounded cancellation polls. Discard output after an error.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_welch_bin(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        recorded_order: &[SourceId],
        budget: HarmonicWorkBudget,
        output: &mut [HarmonicFieldStatistics],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        check_cancel(&mut cancelled)?;
        self.shape(output.len())?;
        let view = source_bin(spectrum, bin, self.omega, &self.sources, recorded_order)?;
        // Includes conversion of the two Cartesian coefficient components.
        self.work(view.segment_count(), 6, false, budget)?;
        for (receiver, out) in output.iter_mut().enumerate() {
            check_cancel(&mut cancelled)?;
            let mut statistics = StatisticsSum::default();
            for factor in view.factors() {
                check_cancel(&mut cancelled)?;
                let mut sums = [Sum::default(); 4];
                for (source, &coefficient) in factor.iter().enumerate() {
                    if source % 256 == 0 {
                        check_cancel(&mut cancelled)?;
                    }
                    let coefficient = peak(view, coefficient)?;
                    let transfer = self.transfers[source * self.receivers.len() + receiver];
                    for component in 0..4 {
                        sums[component].add(transfer[component] * coefficient);
                    }
                }
                let field = sums.map(|sum| sum.finish().scale(core::f64::consts::FRAC_1_SQRT_2));
                statistics.add(field);
            }
            *out = statistics.finish()?;
        }
        check_cancel(&mut cancelled)
    }
}

impl PreparedMicrophoneArray<'_> {
    /// Measured-source covariance of the existing directional microphone bank.
    ///
    /// Returns the full receiver-major Hermitian C_ij [Pa²] IN THIS BIN, with
    /// measured source phase and the original pressure/velocity projections.
    /// A beamformer y=w* m has power w* C w, including cross-receiver terms.
    /// Multiplying this matrix by df again would apply the bandwidth twice.
    /// Source order, physical units and spectral limits follow
    /// `PreparedHarmonicObservation::evaluate_welch_bin`.
    ///
    /// Allocation-free O(segments * (sources * receivers + receivers²)). Each
    /// measured source coefficient is converted once per segment and each
    /// receiver amplitude once per segment, never once per receiver pair.
    pub fn evaluate_welch_bin(
        &mut self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        recorded_order: &[SourceId],
        budget: HarmonicWorkBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&[C64], TransferError> {
        check_cancel(&mut cancelled)?;
        let view = source_bin(
            spectrum,
            bin,
            self.field.omega,
            &self.field.sources,
            recorded_order,
        )?;
        self.field.work(view.segment_count(), 3, true, budget)?;
        self.sums.fill(Sum::default());
        let receivers = self.field.receivers.len();
        for factor in view.factors() {
            check_cancel(&mut cancelled)?;
            let mut sums = [Sum::default(); MAX_RECEIVERS];
            for (source, &coefficient) in factor.iter().enumerate() {
                if source % 256 == 0 {
                    check_cancel(&mut cancelled)?;
                }
                let coefficient = peak(view, coefficient)?;
                for (receiver, sum) in sums[..receivers].iter_mut().enumerate() {
                    sum.add(self.transfers[source * receivers + receiver] * coefficient);
                }
            }
            let mut amplitudes = [C64::ZERO; MAX_RECEIVERS];
            for (amplitude, sum) in amplitudes[..receivers].iter_mut().zip(sums) {
                *amplitude = sum.finish().scale(core::f64::consts::FRAC_1_SQRT_2);
            }
            self.add_mode(&amplitudes[..receivers]);
        }
        check_cancel(&mut cancelled)?;
        self.finish()
    }
}

fn check_cancel(cancelled: &mut impl FnMut() -> bool) -> Result<(), TransferError> {
    if cancelled() {
        Err(TransferError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
