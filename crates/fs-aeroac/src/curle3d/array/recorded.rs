//! Joint physical volume/force recordings and their complete acoustic field.
//!
//! All source-major [Q,Fx,Fy,Fz] channels enter ONE original Welch stream.
//! Cross-source, cross-axis and volume/force interference are retained before
//! pressure/velocity moments are formed. The immutable source layout, physical
//! clock and acoustic reference accompany the factors into propagation.
//!
//! This is spectral observation of accepted physical histories, not a source
//! model, time resampler or waveform synthesizer. Finite-record leakage,
//! unresolved source bandwidth, compact-cell spatial error and the original
//! propagator's domain remain separate. Overlap establishes no independence.

use super::{PreparedCompactSourceObserver, invalid, plane};
use crate::AeroacError;
use crate::curle3d::{CoherenceBudget, CompactFieldStatistics};
use crate::material::MaterialAeroacError;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchStream};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::relaxation::RelaxationState;
use fs_math::c64::C64;

/// The same physical timestamp contract used by the original stress recorder.
/// No second clock or conversion of an adaptive history is introduced.
pub use crate::curle3d::quadrupole::recorded::QuadrupoleSampleClock as SourceSampleClock;

/// Joint bounds on retained source geometry and the original Welch workspace.
#[derive(Clone, Copy, Debug)]
pub struct CompactSourceRecordingBudget {
    pub maximum_sources: usize,
    /// Includes geometry and all original transform/factor buffers. Caller
    /// source fields and allocator metadata are excluded. Four SI channels
    /// are retained per source, including explicitly zero mechanisms.
    pub welch: WelchBudget,
}

/// A bounded recording of prescribed stationary volume and force sources.
///
/// Input channels are real instantaneous Q [m^3/s] and force ON the fluid [N],
/// never phasors or force per span. Every sample retains one original clock,
/// geometry and source order. The numerical acoustic reference is immutable;
/// a molecular propagation law, if requested, is supplied explicitly later.
///
/// Accepted samples allocate nothing. Storage is segment-sized input buffers
/// plus the native spectral factors, not the complete source history or a
/// dense source covariance. Rejected samples do not advance accepted state,
/// including cancellation inside the original transform.
pub struct CompactSourceWelchStream {
    acoustic: FluidAcousticState,
    positions_m: Vec<[f64; 3]>,
    clock: SourceSampleClock,
    stream: WelchStream,
}

impl CompactSourceWelchStream {
    pub fn new(
        acoustic: FluidAcousticState,
        positions_m: &[[f64; 3]],
        clock: SourceSampleClock,
        expected_samples: usize,
        config: WelchConfig,
        budget: CompactSourceRecordingBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        plane::poll(&mut cancelled)?;
        let sources = positions_m.len();
        if sources == 0 || sources > budget.maximum_sources || expected_samples == 0 {
            return Err(invalid(
                "compact recording needs bounded sources and samples",
            ));
        }
        clock.validate_last_index(expected_samples - 1)?;
        let channels = sources
            .checked_mul(4)
            .ok_or_else(|| invalid("compact recording channel count overflow"))?;
        let geometry_bytes = sources
            .checked_mul(size_of::<[f64; 3]>())
            .ok_or_else(|| invalid("compact recording geometry size overflow"))?;
        let mut transform_budget = budget.welch;
        transform_budget.maximum_workspace_bytes = transform_budget
            .maximum_workspace_bytes
            .checked_sub(geometry_bytes)
            .ok_or_else(|| invalid("compact geometry exceeds recording memory budget"))?;
        for (index, position) in positions_m.iter().enumerate() {
            if index.is_multiple_of(256) {
                plane::poll(&mut cancelled)?;
            }
            if position.iter().any(|value| !value.is_finite()) {
                return Err(invalid("compact recording source position is nonfinite"));
            }
        }
        let stream = WelchStream::new(
            channels,
            expected_samples,
            clock.timestep_s(),
            config,
            transform_budget,
            &mut cancelled,
        )
        .map_err(AeroacError::Welch)?;
        let mut positions = Vec::new();
        positions
            .try_reserve_exact(sources)
            .map_err(|_| invalid("compact recording geometry allocation"))?;
        for (index, &position) in positions_m.iter().enumerate() {
            if index.is_multiple_of(256) {
                plane::poll(&mut cancelled)?;
            }
            positions.push(position);
        }
        plane::poll(&mut cancelled)?;
        Ok(Self {
            acoustic,
            positions_m: positions,
            clock,
            stream,
        })
    }

    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    pub const fn clock(&self) -> SourceSampleClock {
        self.clock
    }

    pub fn source_positions_m(&self) -> &[[f64; 3]] {
        &self.positions_m
    }

    pub fn accepted_sample_count(&self) -> usize {
        self.stream.accepted_sample_count()
    }

    pub fn remaining_sample_count(&self) -> usize {
        self.stream.remaining_sample_count()
    }

    pub fn completed_segment_count(&self) -> usize {
        self.stream.completed_segment_count()
    }

    /// Append one complete physical sample without advancing its source solver.
    /// No channels are independently windowed, differentiated or normalized.
    /// Timestamp admission is not a source-bandwidth or aliasing certificate.
    pub fn push_sample(
        &mut self,
        positions_m: &[[f64; 3]],
        sample_time_s: f64,
        volume_force: &[[f64; 4]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        plane::poll(&mut cancelled)?;
        if positions_m.len() != self.positions_m.len()
            || volume_force.len() != self.positions_m.len()
            || self.remaining_sample_count() == 0
        {
            return Err(invalid("compact recording sample shape or declared length"));
        }
        self.clock
            .admit(self.accepted_sample_count(), sample_time_s)?;
        for (index, (actual, expected)) in positions_m.iter().zip(&self.positions_m).enumerate() {
            if index.is_multiple_of(256) {
                plane::poll(&mut cancelled)?;
            }
            if actual != expected {
                return Err(invalid(
                    "compact recording source geometry or order changed",
                ));
            }
        }
        self.stream
            .push_with(
                self.positions_m.len() * 4,
                1,
                |channel, _| volume_force[channel / 4][channel % 4],
                &mut cancelled,
            )
            .map_err(AeroacError::Welch)
    }

    /// Publish the declared complete record without copying spectral factors.
    /// Check remaining_sample_count before consuming an unfinished record.
    pub fn finish(self) -> Result<RecordedCompactSourceSpectrum, AeroacError> {
        Ok(RecordedCompactSourceSpectrum {
            acoustic: self.acoustic,
            positions_m: self.positions_m,
            clock: self.clock,
            spectrum: self.stream.finish().map_err(AeroacError::Welch)?,
        })
    }
}

/// Geometry- and reference-bound joint source factors on their original clock.
///
/// Welch bins carry second moments, not deterministic phases between different
/// frequencies. The retained reference does not independently validate a source
/// model, infer source coherence beyond the record, or identify a molecular law.
pub struct RecordedCompactSourceSpectrum {
    acoustic: FluidAcousticState,
    positions_m: Vec<[f64; 3]>,
    clock: SourceSampleClock,
    spectrum: WelchCrossSpectrum,
}

impl RecordedCompactSourceSpectrum {
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    pub const fn clock(&self) -> SourceSampleClock {
        self.clock
    }

    pub fn source_positions_m(&self) -> &[[f64; 3]] {
        &self.positions_m
    }

    pub fn spectrum(&self) -> &WelchCrossSpectrum {
        &self.spectrum
    }

    /// Prepare the original lossless vector-field transfer at an interior bin.
    /// Existing ideal-wall/wedge methods may subsequently augment this observer;
    /// reflected paths retain their original source's measured coherence.
    pub fn prepare_observer(
        &self,
        observer_m: [f64; 3],
        bin: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedCompactSourceObserver, AeroacError> {
        plane::poll(&mut cancelled)?;
        let frequency_hz = self.bin_frequency(bin)?;
        let k =
            C64::from_re(core::f64::consts::TAU * frequency_hz / self.acoustic.sound_speed_m_s());
        PreparedCompactSourceObserver::prepare_cancellable(
            self.acoustic,
            frequency_hz,
            observer_m,
            &self.positions_m,
            self.positions_m.len(),
            k,
            None,
            &mut cancelled,
        )
    }

    /// Use one explicitly supplied original molecular bulk-response law.
    /// Its equilibrium acoustic reference must match the recording. Original
    /// reservoirs supply absorption and dispersion for pressure AND velocity;
    /// no independently authored attenuation or plane-wave intensity is added.
    pub fn prepare_molecular_observer(
        &self,
        source: RelaxationState,
        observer_m: [f64; 3],
        bin: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedCompactSourceObserver, MaterialAeroacError> {
        plane::poll(&mut cancelled)?;
        if source
            .gas()
            .acoustic_state()
            .map_err(MaterialAeroacError::Material)?
            != self.acoustic
        {
            return Err(invalid("recorded source and molecular acoustic references differ").into());
        }
        let frequency_hz = self.bin_frequency(bin)?;
        let k = source
            .relaxing_wavenumber(core::f64::consts::TAU * frequency_hz)
            .map_err(MaterialAeroacError::Material)?;
        Ok(PreparedCompactSourceObserver::prepare_cancellable(
            self.acoustic,
            frequency_hz,
            observer_m,
            &self.positions_m,
            self.positions_m.len(),
            k,
            Some(source),
            &mut cancelled,
        )?)
    }

    /// Bin-integrated pressure/velocity moments and signed active/reactive
    /// intensity, formed after coherent propagation of every recorded factor.
    /// The original observer enforces exact bin frequency and bounded work.
    pub fn evaluate_bin(
        &self,
        observer: &PreparedCompactSourceObserver,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        plane::poll(&mut cancelled)?;
        if observer.acoustic() != self.acoustic || observer.source_count() != self.positions_m.len()
        {
            return Err(invalid(
                "recorded compact observer changed the medium or source count",
            ));
        }
        for (index, &position) in self.positions_m.iter().enumerate() {
            if index.is_multiple_of(256) {
                plane::poll(&mut cancelled)?;
            }
            if observer.source_position_m(index) != Some(position) {
                return Err(invalid(
                    "recorded compact observer changed source geometry or order",
                ));
            }
        }
        observer.evaluate_welch_bin(&self.spectrum, bin, budget, cancelled)
    }

    fn bin_frequency(&self, bin: usize) -> Result<f64, AeroacError> {
        let frequency = self
            .spectrum
            .frequency_hz(bin)
            .map_err(AeroacError::Welch)?;
        self.spectrum
            .peak_bin_hz(bin, frequency)
            .map_err(AeroacError::Welch)?;
        Ok(frequency)
    }
}

#[cfg(test)]
mod tests;
