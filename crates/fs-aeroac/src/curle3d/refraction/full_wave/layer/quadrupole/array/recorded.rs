//! Original acquired stress histories through their declared layer observation.

use super::*;
use crate::curle3d::quadrupole::recorded::RecordedQuadrupoleSpectrum;

impl RecordedQuadrupoleSpectrum {
    /// Bind an original interior Welch bin and original source geometry to a
    /// transmitted, reflected, coherent total or interior layer observation.
    /// The recording owns its clock, stress units, window and sample history;
    /// the layer owns the three admitted original fluid states and faces.
    pub fn prepare_layered_observer(
        &self,
        layer: PlanarFluidLayer,
        field_kind: LayerQuadrupoleField,
        observer_m: [f64; 3],
        bin: usize,
        limits: FullWaveQuadrupoleSourceLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedLayeredFullWaveQuadrupoleObserver, AeroacError> {
        poll(&mut cancelled)?;
        let frequency_hz = self
            .spectrum()
            .frequency_hz(bin)
            .map_err(AeroacError::Welch)?;
        self.spectrum()
            .peak_bin_hz(bin, frequency_hz)
            .map_err(AeroacError::Welch)?;
        PreparedLayeredFullWaveQuadrupoleObserver::new(
            layer,
            field_kind,
            frequency_hz,
            observer_m,
            self.source_positions_m(),
            limits,
            &mut cancelled,
        )
    }

    /// Preserve source ordering and positions before applying the original
    /// joint Welch factors. sqrt(2 df) is applied once by the native peak-bin
    /// owner. No new density, dV, phase origin, ray delay, FFT, covariance matrix,
    /// or inferred independence is introduced. Output is bin-integrated power.
    pub fn evaluate_layered_bin(
        &self,
        observer: &PreparedLayeredFullWaveQuadrupoleObserver,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        poll(&mut cancelled)?;
        observer
            .view()
            .admit_count(self.source_positions_m().len())?;
        for (index, &position) in self.source_positions_m().iter().enumerate() {
            if index.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            if observer.source_position_m(index) != Some(position) {
                return Err(invalid(
                    "layered observer changed recorded stress geometry or order",
                ));
            }
        }
        observer.evaluate_welch_bin(self.spectrum(), bin, budget, cancelled)
    }
}
