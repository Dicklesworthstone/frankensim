//! Original acquired inside-layer stress histories and their exterior radiation.

use super::*;
use crate::curle3d::quadrupole::recorded::RecordedQuadrupoleSpectrum;

impl RecordedQuadrupoleSpectrum {
    /// Bind one original interior Welch bin and original source positions to an
    /// inside-layer source observation in the incident exterior. The recording
    /// owns the source clock, physical stress moments, window and sample history;
    /// the native layer owner admits the actual source and receiving fluids.
    pub fn prepare_interior_source_layered_observer(
        &self,
        layer: PlanarFluidLayer,
        observer_m: [f64; 3],
        bin: usize,
        limits: FullWaveQuadrupoleSourceLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedInteriorSourceFullWaveQuadrupoleObserver, AeroacError> {
        poll(&mut cancelled)?;
        let frequency_hz = self
            .spectrum()
            .frequency_hz(bin)
            .map_err(AeroacError::Welch)?;
        self.spectrum()
            .peak_bin_hz(bin, frequency_hz)
            .map_err(AeroacError::Welch)?;
        PreparedInteriorSourceFullWaveQuadrupoleObserver::new(
            layer,
            frequency_hz,
            observer_m,
            self.source_positions_m(),
            limits,
            &mut cancelled,
        )
    }

    /// Preserve original source order and geometry before applying the joint
    /// factors. Native sqrt(2 df) normalization occurs once; output is bin power.
    /// No additional density, cell volume, phase origin or ray delay is applied,
    /// and tensor/cell coherence is never replaced by inferred independence.
    pub fn evaluate_interior_source_layered_bin(
        &self,
        observer: &PreparedInteriorSourceFullWaveQuadrupoleObserver,
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
                    "inside-layer observer changed recorded stress geometry or order",
                ));
            }
        }
        observer.evaluate_welch_bin(self.spectrum(), bin, budget, cancelled)
    }
}
