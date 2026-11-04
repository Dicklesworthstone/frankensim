//! Recorded source coherence through the original loaded network operator.

use super::*;
use fs_fft::welch::{WelchCrossSpectrum, WelchPeakBin};

impl PreparedNetworkSources {
    /// Propagate one measured interior frequency bin through the original
    /// pressure, flow, internal-source and dissipative-source transfer banks.
    /// All cross-port p/p, U/U and p/U moments and signed work summaries use
    /// the same contractions as [`Self::coherence_into`]. No graph solve,
    /// constitutive query, dense source covariance or factor copy is repeated.
    ///
    /// Recorded order must exactly match source labels AND physical quantities:
    /// prescribed volume flow is m³/s; open-circuit source pressure is Pa before
    /// its entire original impedance. Source resistance and external loading
    /// are already contained in this preparation. Their mutual interference
    /// and signed energy exchange are retained, without clipping or repair.
    ///
    /// The shared bin view applies sqrt(2*df) once, followed by the original
    /// peak time-average half once. Results are bin moments and watts, not
    /// densities per Hz. DC/Nyquist and mismatched physical frequencies refuse.
    /// Segment factors may overlap; this does not infer independence,
    /// cross-frequency coherence, statistical confidence or stationarity.
    ///
    /// Allocation-free; at most 4096 recorded factors are admitted under the
    /// existing product ceiling. Cancellation polls each factor, with at most
    /// the original bounded 32-port contraction between polls. Any refusal
    /// invalidates output and its getters, including after a previous success.
    #[allow(clippy::too_many_arguments)]
    pub fn observe_welch_bin_into(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        recorded_order: &[NetworkSourceCoordinate],
        budget: NetworkCoherenceBudget,
        output: &mut NetworkSourceStatistics,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), DuctError> {
        output.valid = false;
        check_cancel(&mut cancelled, 0)?;
        let view = self.recorded_bin(spectrum, bin, recorded_order)?;
        self.admit_statistics(output, view.segment_count(), false, budget, true)?;
        let count = self.coordinates.len();
        let mut excitation = [C64::ZERO; MAX_NETWORK_PORTS];
        output.clear();
        for (segment, factor) in view.factors().enumerate() {
            check_cancel(&mut cancelled, segment)?;
            for (out, &value) in excitation[..count].iter_mut().zip(factor) {
                let value = view
                    .peak_factor(value)
                    .map_err(|_| bad("recorded network source peak is unrepresentable"))?;
                *out = C64::new(value.re, value.im);
            }
            let sample = self.sample(&excitation[..count], None)?;
            self.accumulate(&sample, &excitation[..count], output)?;
        }
        check_cancel(&mut cancelled, view.segment_count())?;
        self.finish_statistics(output)
    }

    pub(in crate::network) fn recorded_bin<'a>(
        &self,
        spectrum: &'a WelchCrossSpectrum,
        bin: usize,
        recorded_order: &[NetworkSourceCoordinate],
    ) -> Result<WelchPeakBin<'a>, DuctError> {
        let view = spectrum
            .peak_bin_rad_s(bin, self.omega)
            .map_err(|_| bad("recorded network bin must match the positive physical frequency"))?;
        if recorded_order != self.coordinates || view.channel_count() != self.coordinates.len() {
            return Err(bad(
                "recorded network labels, SI source quantities or order differ",
            ));
        }
        Ok(view)
    }
}

fn check_cancel(
    cancelled: &mut impl FnMut() -> bool,
    completed_steps: usize,
) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps })
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
