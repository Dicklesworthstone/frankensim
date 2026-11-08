//! Original joint Q/F/stress factors through the native wall transfers.

use super::super::PreparedImpedanceWallQuadrupoleObserver;
use super::super::recorded::{same_geometry, scope};
use super::*;
use crate::curle3d::statistics::{ComplexSum, StatisticsSum};
use fs_fft::welch::{WelchCrossSpectrum, WelchPeakBin};

impl PreparedImpedanceWallCompactSourceObserver {
    /// Observe one interior bin of original synchronized `[Q,Fx,Fy,Fz]`
    /// source blocks, in the prepared position order. Q is injected volume
    /// flow [m³/s]; F is force ON fluid [N]. Missing mechanisms are zero
    /// channels. The caller binds the original recording geometry and order.
    ///
    /// Each original Welch factor contributes one coherent pressure/vector-
    /// velocity sum before second moments are formed. The existing bin owner
    /// applies sqrt(2 df) exactly once. Results are bin-integrated moments,
    /// retaining the original window, detrending, leakage and sampling limits.
    /// No independent-source assumption, covariance expansion, extra delay,
    /// FFT or propagation quadrature is introduced. Application is allocation-
    /// free, O(sources * segments), with bounded work and cancellation.
    pub fn evaluate_welch_bin(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        let mut statistics = StatisticsSum::default();
        self.accumulate_welch_bin(None, spectrum, bin, budget, &mut statistics, &mut cancelled)?;
        statistics.finish(scope(self.field_kind()))
    }

    /// Observe a single original joint record of compact AND volume-stress
    /// sources: 4N source-major `[Q,Fx,Fy,Fz]` channels followed by 6M physical
    /// `[xx,yy,zz,xy,xz,yz]` stress moments [N m]. Both banks retain their own
    /// prepared source order. Stress shear is not pre-doubled; density and
    /// cell volume are already included in each moment.
    ///
    /// The banks must share the same fluid, wall geometry AND impedance,
    /// receiver, exact bin frequency and reflected/total field. Both groups
    /// enter one DD pressure and velocity sum within every original factor.
    /// Cross-source, Q/F/stress and direct/reflected interference survives,
    /// including nonzero velocity at a pressure null. Separate spectra cannot
    /// supply missing cross-mechanism phase. Physical mechanisms must be
    /// disjoint; do not recount a volume represented by a permeable surface.
    ///
    /// The combined source/factor budget is checked before propagation.
    /// These are bin-integrated statistics, not deterministic cross-frequency
    /// phases, a time waveform, a confidence interval or inferred SGS sound.
    pub fn evaluate_welch_bin_with_quadrupoles(
        &self,
        quadrupoles: &PreparedImpedanceWallQuadrupoleObserver,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        let mut statistics = StatisticsSum::default();
        self.accumulate_welch_bin(
            Some(quadrupoles),
            spectrum,
            bin,
            budget,
            &mut statistics,
            &mut cancelled,
        )?;
        statistics.finish(scope(self.field_kind()))
    }

    /// The band owner reuses this same factor reduction, accumulating distinct
    /// resolved bins into scratch statistics without rounding each bin first.
    pub(in crate::curle3d::refraction::full_wave::wall) fn accumulate_welch_bin(
        &self,
        quadrupoles: Option<&PreparedImpedanceWallQuadrupoleObserver>,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        statistics: &mut StatisticsSum,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        poll(cancelled)?;
        if let Some(stress) = quadrupoles {
            if !same_geometry(self.wall(), stress.wall())
                || self.wall().impedance_pa_s_m() != stress.wall().impedance_pa_s_m()
                || self.observer_position_m() != stress.observer_position_m()
                || self.frequency_hz().to_bits() != stress.frequency_hz().to_bits()
                || self.field_kind() != stress.field_kind()
            {
                return Err(invalid(
                    "recorded wall source banks need one fluid, geometry, impedance, receiver, frequency and field",
                ));
            }
        }
        let stress_sources = quadrupoles.map_or(0, |q| q.source_count());
        let compact_channels = self
            .source_count()
            .checked_mul(4)
            .ok_or_else(|| invalid("recorded wall compact channel count overflow"))?;
        let channels = stress_sources
            .checked_mul(6)
            .and_then(|n| n.checked_add(compact_channels))
            .ok_or_else(|| invalid("recorded wall mixed channel count overflow"))?;
        let sources = self
            .source_count()
            .checked_add(stress_sources)
            .ok_or_else(|| invalid("recorded wall mixed source count overflow"))?;
        let measured = spectrum
            .peak_bin_hz(bin, self.frequency_hz())
            .map_err(AeroacError::Welch)?;
        if channels != measured.channel_count()
            || measured.segment_count() > budget.max_modes
            || sources
                .checked_mul(measured.segment_count())
                .is_none_or(|work| work > budget.max_source_mode_pairs)
        {
            return Err(invalid(
                "recorded wall source shape or combined coherence budget",
            ));
        }
        for factor in measured.factors() {
            poll(cancelled)?;
            let (compact, stress) = factor.split_at(compact_channels);
            let mut sums = [ComplexSum::default(); 4];
            self.accumulate_welch_factor(compact, measured, &mut sums, cancelled)?;
            if let Some(quadrupoles) = quadrupoles {
                quadrupoles.accumulate_welch_factor(stress, measured, &mut sums, cancelled)?;
            }
            statistics.add(sums.map(ComplexSum::value))?;
        }
        poll(cancelled)
    }

    fn accumulate_welch_factor(
        &self,
        factor: &[fs_fft::C64],
        measured: WelchPeakBin<'_>,
        sums: &mut [ComplexSum; 4],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        for (index, (transfer, channels)) in self
            .transfers
            .iter()
            .zip(factor.chunks_exact(4))
            .enumerate()
        {
            if index.is_multiple_of(256) {
                poll(cancelled)?;
            }
            let mut values = [C64::ZERO; 4];
            for (coefficient, value) in channels.iter().zip(&mut values) {
                let peak = measured
                    .peak_factor(*coefficient)
                    .map_err(AeroacError::Welch)?;
                *value = C64::new(peak.re, peak.im);
            }
            transfer.accumulate(
                CompactSourcePhasor {
                    volume_flow_m3_s: values[0],
                    force_on_fluid_n: [values[1], values[2], values[3]],
                },
                sums,
            )?;
        }
        Ok(())
    }
}
