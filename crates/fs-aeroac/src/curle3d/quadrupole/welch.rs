//! Recorded volume stresses through the original coherent quadrupole transfer.

use super::*;
use fs_fft::welch::WelchCrossSpectrum;

impl PreparedQuadrupoleObserver {
    /// Propagate one interior bin of synchronized physical stress histories.
    ///
    /// Channels are source-major `[Qxx, Qyy, Qzz, Qxy, Qxz, Qyz]` in N m,
    /// in the prepared positions' exact order and Cartesian frame. Shear
    /// entries are physical tensor entries, not pre-doubled. Record the
    /// instantaneous stress before estimating its spectrum: a spectrum of
    /// velocity is not a spectrum of the quadratic Reynolds-stress source.
    ///
    /// Each original Welch factor is propagated coherently across every cell
    /// and tensor component before pressure/velocity second moments are
    /// formed. Relative phases and cross-component/cross-cell terms survive.
    /// The original one-sided /Hz factors receive `sqrt(2 df)` exactly once;
    /// outputs are bin-integrated moments, not PSD/Hz or coherent-gain tone
    /// amplitudes. No extra density, cell volume or propagation phase is added.
    ///
    /// This retains the estimator's window, detrending, leakage, finite-record
    /// and source-model limits. Overlapping segments are not claimed to be
    /// independent. Summing distinct bins requires a stationary resolved
    /// spectrum; this supplies no cross-frequency beats, deterministic signal
    /// reconstruction, confidence interval or unresolved SGS completion.
    ///
    /// Allocation-free O(cells * segments), using prepared transfer columns
    /// and deterministic double-double reductions. There is no FFT, dense
    /// source covariance or transcendental work in this observation. Source
    /// count, frequency and work admission precede propagation; cancellation
    /// is checked at least every 256 cells and between recorded factors.
    pub fn evaluate_welch_bin(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        poll(&mut cancelled)?;
        let measured = spectrum
            .peak_bin_hz(bin, self.frequency_hz)
            .map_err(AeroacError::Welch)?;
        if self.source_count().checked_mul(6) != Some(measured.channel_count()) {
            return Err(invalid(
                "recorded quadrupoles require six stress channels per prepared source",
            ));
        }
        if measured.segment_count() > budget.max_modes
            || self
                .source_count()
                .checked_mul(measured.segment_count())
                .is_none_or(|work| work > budget.max_source_mode_pairs)
        {
            return Err(invalid("recorded quadrupole coherence exceeds work budget"));
        }
        let mut statistics = StatisticsSum::default();
        for factor in measured.factors() {
            poll(&mut cancelled)?;
            let mut sums = [ComplexSum::default(); 4];
            for (index, (transfer, channels)) in self
                .transfers
                .iter()
                .zip(factor.chunks_exact(6))
                .enumerate()
            {
                if index.is_multiple_of(256) {
                    poll(&mut cancelled)?;
                }
                for (entry, coefficient) in channels.iter().enumerate() {
                    let peak = measured
                        .peak_factor(*coefficient)
                        .map_err(AeroacError::Welch)?;
                    let value = C64::new(peak.re, peak.im);
                    if value != C64::ZERO {
                        for (row, sum) in transfer.field.iter().zip(&mut sums) {
                            sum.add(row[entry] * value)?;
                        }
                    }
                }
            }
            statistics.add(sums.map(ComplexSum::value))?;
        }
        poll(&mut cancelled)?;
        statistics.finish(QUADRUPOLE_SCOPE)
    }
}
