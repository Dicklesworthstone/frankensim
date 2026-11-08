//! Original recorded volume-stress cells beside a passive impedance wall.

use super::*;
use crate::curle3d::quadrupole::recorded::RecordedQuadrupoleSpectrum;
use fs_fft::welch::WelchCrossSpectrum;

impl PreparedImpedanceWallQuadrupoleObserver {
    /// Propagate one original interior Welch bin through the cached wall field.
    ///
    /// Source-major channels are physical [xx,yy,zz,xy,xz,yz] stress moments
    /// in N m, in the prepared cell order. Every original factor keeps its cell,
    /// tensor and direct/reflected phases before pressure/velocity moments are
    /// formed. The shared peak-bin owner applies sqrt(2 df) exactly once.
    ///
    /// The recorded frequency must exactly match this prepared wall response.
    /// Results are bin-integrated second moments, retaining the original window,
    /// finite-record and aliasing limitations; they supply no cross-frequency
    /// phase, waveform continuation or independent-cell assumption.
    ///
    /// Application is allocation-free O(cells * segments), with explicit work
    /// limits and bounded cancellation. It reuses the original transfer and
    /// factor owners without another FFT, quadrature, source volume or delay.
    /// Use the recording-aware method to bind original positions and ordering.
    pub fn evaluate_welch_bin(
        &self,
        spectrum: &WelchCrossSpectrum,
        bin: usize,
        budget: CoherenceBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view()
            .evaluate_welch_bin(spectrum, bin, budget, cancelled)
    }
}

impl RecordedQuadrupoleSpectrum {
    /// Bind the original recorded cells and interior-bin frequency to their
    /// reflected-only pressure and vector velocity at a passive wall. Sources
    /// lie strictly in the fluid; the receiver may be on the wall or coincide
    /// with a cell for its finite reflected self-field.
    ///
    /// The supplied impedance must be evaluated by its original material owner
    /// at this recorded bin's frequency. Preparation does not infer a broadband
    /// material law or replace finite impedance with a constant image strength.
    pub fn prepare_impedance_wall_reflected_observer(
        &self,
        wall: PlanarImpedanceWall,
        observer_m: [f64; 3],
        bin: usize,
        limits: FullWaveQuadrupoleSourceLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedImpedanceWallQuadrupoleObserver, AeroacError> {
        self.prepare_impedance_wall_observer(
            wall,
            observer_m,
            bin,
            limits,
            ImpedanceWallFieldKind::Reflected,
            cancelled,
        )
    }

    /// Bind coherent direct-plus-reflected observation to the same original
    /// stress record. Every cell's paths interfere inside each measured factor
    /// before statistical moments are formed. Singular direct self-fields refuse.
    /// The caller supplies the wall impedance at the original bin frequency.
    pub fn prepare_impedance_wall_total_observer(
        &self,
        wall: PlanarImpedanceWall,
        observer_m: [f64; 3],
        bin: usize,
        limits: FullWaveQuadrupoleSourceLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedImpedanceWallQuadrupoleObserver, AeroacError> {
        self.prepare_impedance_wall_observer(
            wall,
            observer_m,
            bin,
            limits,
            ImpedanceWallFieldKind::Total,
            cancelled,
        )
    }

    fn prepare_impedance_wall_observer(
        &self,
        wall: PlanarImpedanceWall,
        observer_m: [f64; 3],
        bin: usize,
        limits: FullWaveQuadrupoleSourceLimits,
        field_kind: ImpedanceWallFieldKind,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedImpedanceWallQuadrupoleObserver, AeroacError> {
        poll(&mut cancelled)?;
        let frequency_hz = self
            .spectrum()
            .frequency_hz(bin)
            .map_err(AeroacError::Welch)?;
        self.spectrum()
            .peak_bin_hz(bin, frequency_hz)
            .map_err(AeroacError::Welch)?;
        PreparedImpedanceWallQuadrupoleObserver::prepare(
            wall,
            field_kind,
            frequency_hz,
            observer_m,
            self.source_positions_m(),
            limits,
            &mut cancelled,
        )
    }

    /// Apply the original joint factors only to their original geometry and
    /// ordering, for reflected-only or coherent total wall observations.
    ///
    /// Equal cell counts do not establish identity: relocated or reordered
    /// sources refuse before reduction. The shared recorded-bin owner then
    /// enforces exact frequency, six channels per cell and source-factor work.
    /// Original stress units and sample clock remain unchanged.
    pub fn evaluate_impedance_wall_bin(
        &self,
        observer: &PreparedImpedanceWallQuadrupoleObserver,
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
                    "impedance wall observer changed recorded stress geometry or order",
                ));
            }
        }
        observer.evaluate_welch_bin(self.spectrum(), bin, budget, cancelled)
    }
}
