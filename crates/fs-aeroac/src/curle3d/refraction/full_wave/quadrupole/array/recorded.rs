//! The original joint volume-stress spectrum through refracting transfers.

use super::*;
use crate::curle3d::quadrupole::recorded::RecordedQuadrupoleSpectrum;
use fs_fft::welch::WelchCrossSpectrum;

impl PreparedFullWaveQuadrupoleObserver {
    /// Propagate one recorded interior Welch bin without losing cell coherence.
    ///
    /// Source-major channels are the original [xx,yy,zz,xy,xz,yz] stress moments
    /// in N m. The native peak-bin owner applies sqrt(2 df) exactly once to its
    /// one-sided /Hz factors. The output is a bin-integrated second moment,
    /// retaining the original window/detrending/finite-record limitations.
    /// This does not infer independent cells, stationary cross-frequency phase,
    /// source calibration, a waveform, extra cell volume or another ray delay.
    ///
    /// Application is allocation-free O(cells*segments), with bounded polling.
    /// Geometry identity is supplied by the recording-aware method below.
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
                "refracted recorded stresses require six channels per cell",
            ));
        }
        self.admit_modes(measured.segment_count(), budget)?;
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
        statistics.finish(self.scope())
    }
}

impl RecordedQuadrupoleSpectrum {
    /// Prepare refracting propagation at one original recorded frequency and
    /// source layout. The existing recording owns its clock, stress normalization
    /// and positions; the supplied interface owns both actual fluid states.
    pub fn prepare_refracted_observer(
        &self,
        interface: PlanarFluidInterface,
        observer_m: [f64; 3],
        bin: usize,
        limits: FullWaveQuadrupoleSourceLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedFullWaveQuadrupoleObserver, AeroacError> {
        poll(&mut cancelled)?;
        let frequency_hz = self
            .spectrum()
            .frequency_hz(bin)
            .map_err(AeroacError::Welch)?;
        self.spectrum()
            .peak_bin_hz(bin, frequency_hz)
            .map_err(AeroacError::Welch)?;
        PreparedFullWaveQuadrupoleObserver::new(
            interface,
            frequency_hz,
            observer_m,
            self.source_positions_m(),
            limits,
            &mut cancelled,
        )
    }

    /// Propagate the same recorded factors only through their original geometry.
    /// Same-sized but reordered or relocated cells refuse before any reduction.
    pub fn evaluate_refracted_bin(
        &self,
        observer: &PreparedFullWaveQuadrupoleObserver,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        poll(&mut cancelled)?;
        observer.admit_count(self.source_positions_m().len())?;
        for (i, &position) in self.source_positions_m().iter().enumerate() {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            if observer.source_position_m(i) != Some(position) {
                return Err(invalid(
                    "refracted observer changed recorded stress geometry or order",
                ));
            }
        }
        observer.evaluate_welch_bin(self.spectrum(), bin, budget, cancelled)
    }
}
