//! Measured source factors use the original prepared pressure/velocity paths.

use super::*;
use fs_fft::welch::{WelchCrossSpectrum, WelchPeakBin};

impl PreparedCompactSourceObserver {
    /// Propagate one measured interior Welch bin through the original medium.
    ///
    /// Channels are source-major `[Q, Fx, Fy, Fz]`: physical volume injection
    /// in m³/s followed by force ON the fluid in N for each prepared position.
    /// All channels must share the original sampling clock. A missing mechanism
    /// is an explicitly zero channel; no span conversion or force is inferred.
    ///
    /// Each recorded factor is propagated coherently through the same cached
    /// pressure/velocity coefficients, including any molecular dispersion or
    /// ideal-plane image paths. Source, force-axis and volume/force cross terms
    /// survive. Images never become additional independent source channels.
    ///
    /// The one-sided /Hz factors are scaled by sqrt(2 df) exactly once before
    /// the existing peak-phasor statistics. Results are bin-integrated second
    /// moments, not PSD/Hz or coherent-gain tone amplitudes. Window leakage and
    /// the original source/propagation scope remain. Adding bins assumes a
    /// stationary resolved spectrum; this supplies no cross-frequency beats,
    /// confidence interval or statistical independence of overlapping segments.
    ///
    /// Allocation-free O(sources * segments), with source-boundary cancellation
    /// and no FFT, Green-function evaluation or dense source covariance.
    pub fn evaluate_welch_bin(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        plane::poll(&mut cancelled)?;
        let measured = admit(
            spectrum,
            bin,
            self.frequency_hz,
            self.source_count(),
            4,
            1,
            budget,
        )?;
        let mut statistics = StatisticsSum::default();
        for factor in measured.factors() {
            plane::poll(&mut cancelled)?;
            let mut sums = [ComplexSum::default(); 4];
            for (transfer, channels) in self.transfers.iter().zip(factor.chunks_exact(4)) {
                plane::poll(&mut cancelled)?;
                transfer.accumulate(source(channels, measured)?, &mut sums)?;
            }
            statistics.add(sums.map(ComplexSum::value))?;
        }
        plane::poll(&mut cancelled)?;
        statistics.finish(self.scope())
    }
}

/// Same measured-bin admission for compact, half-space and convected owners.
/// `paths_per_source` counts separately applied paths, not folded image paths.
pub(super) fn admit(
    spectrum: &WelchCrossSpectrum,
    bin: usize,
    frequency_hz: f64,
    sources: usize,
    channels_per_source: usize,
    paths_per_source: usize,
    budget: CoherenceBudget,
) -> Result<WelchPeakBin<'_>, AeroacError> {
    let measured = spectrum
        .peak_bin_hz(bin, frequency_hz)
        .map_err(|_| invalid("recorded compact source bin/frequency"))?;
    if sources.checked_mul(channels_per_source) != Some(measured.channel_count()) {
        return Err(invalid(
            "Welch source channels must match the prepared compact observer",
        ));
    }
    if measured.segment_count() > budget.max_modes
        || sources
            .checked_mul(measured.segment_count())
            .and_then(|work| work.checked_mul(paths_per_source))
            .is_none_or(|work| work > budget.max_source_mode_pairs)
    {
        return Err(invalid(
            "recorded compact source coherence exceeds work budget",
        ));
    }
    Ok(measured)
}

pub(super) fn source(
    channels: &[fs_fft::C64],
    bin: WelchPeakBin<'_>,
) -> Result<CompactSourcePhasor, AeroacError> {
    let value = |i: usize| -> Result<C64, AeroacError> {
        let peak = bin
            .peak_factor(channels[i])
            .map_err(|_| invalid("recorded compact source peak factor is unrepresentable"))?;
        Ok(C64::new(peak.re, peak.im))
    };
    Ok(CompactSourcePhasor {
        volume_flow_m3_s: value(0)?,
        force_on_fluid_n: [value(1)?, value(2)?, value(3)?],
    })
}
