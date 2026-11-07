//! Clock- and geometry-bound acquisition of physical volume-stress histories.
//!
//! The six instantaneous stress moments of every cell share one original
//! sampling clock. The existing joint Welch owner retains their relative phase
//! and applies its window, detrending and normalization exactly once. Nothing
//! reconstructs stress from a velocity spectrum or invents independent cells.

use super::{CoherenceBudget, CompactFieldStatistics, PreparedQuadrupoleObserver, poll};
use crate::AeroacError;
use crate::material::MaterialAeroacError;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchStream};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::relaxation::RelaxationState;
use fs_math::dd::Dd;

/// Uniform physical source clock, with an explicit timestamp admission error.
///
/// The tolerance admits only the caller's declared timestamp discrepancy. It
/// does not resample a nonuniform history or bound its spectral aliasing error.
/// A truly variable-step solver must first supply a qualified uniform record.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadrupoleSampleClock {
    epoch_s: f64,
    timestep_s: f64,
    maximum_timestamp_error_s: f64,
}

impl QuadrupoleSampleClock {
    /// The epoch is the time of sample zero, not a newly invented phase origin.
    /// Require 0 <= timestamp error < dt/4; zero requests exact timestamps.
    pub fn new(
        epoch_s: f64,
        timestep_s: f64,
        maximum_timestamp_error_s: f64,
    ) -> Result<Self, AeroacError> {
        if !epoch_s.is_finite()
            || !timestep_s.is_finite()
            || timestep_s <= 0.0
            || !maximum_timestamp_error_s.is_finite()
            || maximum_timestamp_error_s < 0.0
            || maximum_timestamp_error_s >= 0.25 * timestep_s
        {
            return Err(invalid("quadrupole source clock or timestamp error"));
        }
        let clock = Self {
            epoch_s,
            timestep_s,
            maximum_timestamp_error_s,
        };
        clock.validate_last_index(1)?;
        Ok(clock)
    }

    pub const fn epoch_s(self) -> f64 {
        self.epoch_s
    }

    pub const fn timestep_s(self) -> f64 {
        self.timestep_s
    }

    pub const fn maximum_timestamp_error_s(self) -> f64 {
        self.maximum_timestamp_error_s
    }

    /// Compute from the original epoch and integer index, without accumulated
    /// addition drift. Recording preparation separately admits the timestamp
    /// resolution of the complete requested range.
    pub fn time_at(self, sample_index: usize) -> Result<f64, AeroacError> {
        // Every admitted integer is exactly representable by f64 on both
        // 32- and 64-bit targets. Do not silently lose a large sample count.
        if (sample_index as u128) > (1u128 << 53) {
            return Err(invalid("quadrupole sample index is not representable"));
        }
        let time = Dd::from_f64(self.epoch_s)
            + Dd::from_f64(self.timestep_s) * Dd::from_f64(sample_index as f64);
        let value = time.to_f64();
        if !time.hi.is_finite() || !time.lo.is_finite() || !value.is_finite() {
            return Err(invalid("quadrupole sample timestamp is not representable"));
        }
        Ok(value)
    }

    pub(super) fn validate_last_index(self, last: usize) -> Result<(), AeroacError> {
        let largest = self.epoch_s.abs().max(self.time_at(last)?.abs());
        let spacing = largest.next_up() - largest;
        // Checking one adjacent pair is insufficient: dt below one ULP can
        // alternate between distinct and coincident rounded timestamps.
        if self.timestep_s < spacing
            || (last != 0 && self.time_at(last)? <= self.time_at(last - 1)?)
        {
            return Err(invalid(
                "quadrupole timestep is unresolved at the source epoch",
            ));
        }
        Ok(())
    }

    pub(super) fn admit(self, index: usize, time_s: f64) -> Result<(), AeroacError> {
        let expected = self.time_at(index)?;
        if !time_s.is_finite()
            || (Dd::from_f64(time_s) - Dd::from_f64(expected))
                .to_f64()
                .abs()
                > self.maximum_timestamp_error_s
        {
            return Err(invalid(
                "quadrupole sample is off the declared physical clock",
            ));
        }
        Ok(())
    }
}

/// Bounds the source geometry and the original streaming transform together.
#[derive(Clone, Copy, Debug)]
pub struct QuadrupoleRecordingBudget {
    pub maximum_sources: usize,
    /// Workspace bytes include retained positions plus all original Welch
    /// buffers/factors. Caller-owned flow fields and allocator metadata are not
    /// counted. Channels are exactly six physical tensor entries per source.
    pub welch: WelchBudget,
}

/// Streaming volume histories with fixed physical source order and geometry.
///
/// Source-major entries are xx, yy, zz, xy, xz, yz in N m. Memory follows the
/// original Welch owner: segment-sized input buffers and retained spectral
/// factors, never the full flow history or a dense source covariance. An
/// accepted sample makes no allocation. Refusal advances neither sample count
/// nor completed segment count, including cancellation during a transform.
pub struct QuadrupoleWelchStream {
    source_positions_m: Vec<[f64; 3]>,
    clock: QuadrupoleSampleClock,
    stream: WelchStream,
}

impl QuadrupoleWelchStream {
    pub fn new(
        source_positions_m: &[[f64; 3]],
        clock: QuadrupoleSampleClock,
        expected_samples: usize,
        config: WelchConfig,
        budget: QuadrupoleRecordingBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        let sources = source_positions_m.len();
        let channels = sources.checked_mul(6);
        let position_bytes = sources.checked_mul(core::mem::size_of::<[f64; 3]>());
        if sources == 0 || sources > budget.maximum_sources || expected_samples == 0 {
            return Err(invalid(
                "quadrupole recording needs bounded sources and samples",
            ));
        }
        clock.validate_last_index(expected_samples - 1)?;
        let channels = channels.ok_or_else(|| invalid("quadrupole channel count overflow"))?;
        let position_bytes =
            position_bytes.ok_or_else(|| invalid("quadrupole geometry storage overflow"))?;
        let mut transform_budget = budget.welch;
        transform_budget.maximum_workspace_bytes = transform_budget
            .maximum_workspace_bytes
            .checked_sub(position_bytes)
            .ok_or_else(|| invalid("quadrupole geometry exceeds recording memory budget"))?;
        // Admission before allocation keeps invalid geometry from consuming a
        // transform plan or factor bank.
        for (i, position) in source_positions_m.iter().enumerate() {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            if position.iter().any(|value| !value.is_finite()) {
                return Err(invalid("quadrupole source position is nonfinite"));
            }
        }
        let stream = WelchStream::new(
            channels,
            expected_samples,
            clock.timestep_s,
            config,
            transform_budget,
            &mut cancelled,
        )
        .map_err(AeroacError::Welch)?;
        let mut positions = Vec::new();
        positions
            .try_reserve_exact(sources)
            .map_err(|_| invalid("quadrupole source-position allocation"))?;
        for (i, &position) in source_positions_m.iter().enumerate() {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            positions.push(position);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            source_positions_m: positions,
            clock,
            stream,
        })
    }

    pub const fn clock(&self) -> QuadrupoleSampleClock {
        self.clock
    }

    pub fn source_positions_m(&self) -> &[[f64; 3]] {
        &self.source_positions_m
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

    /// Observe already accepted physical moments without changing the solver.
    /// Positions must retain the prepared order exactly. The caller supplies
    /// moments from this time; neither a stale flow sample nor adequate source
    /// bandwidth can be detected from its timestamp alone. There is no implicit
    /// subsampling: supply every admitted source-clock sample.
    pub fn push_sample(
        &mut self,
        source_positions_m: &[[f64; 3]],
        sample_time_s: f64,
        stress_moments_n_m: &[[f64; 6]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        poll(&mut cancelled)?;
        if source_positions_m.len() != self.source_positions_m.len()
            || stress_moments_n_m.len() != self.source_positions_m.len()
            || self.remaining_sample_count() == 0
        {
            return Err(invalid("quadrupole sample shape or declared record length"));
        }
        self.clock
            .admit(self.accepted_sample_count(), sample_time_s)?;
        for (i, (actual, expected)) in source_positions_m
            .iter()
            .zip(&self.source_positions_m)
            .enumerate()
        {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            if actual != expected {
                return Err(invalid("quadrupole source geometry or order changed"));
            }
        }
        self.stream
            .push_with(
                self.source_positions_m.len() * 6,
                1,
                |channel, _| stress_moments_n_m[channel / 6][channel % 6],
                &mut cancelled,
            )
            .map_err(AeroacError::Welch)
    }

    /// Publish the complete declared record without copying its factors.
    /// Check remaining_sample_count before consuming an incomplete record.
    pub fn finish(self) -> Result<RecordedQuadrupoleSpectrum, AeroacError> {
        Ok(RecordedQuadrupoleSpectrum {
            source_positions_m: self.source_positions_m,
            clock: self.clock,
            spectrum: self.stream.finish().map_err(AeroacError::Welch)?,
        })
    }
}

/// Original joint stress spectrum with its permanent physical source layout.
///
/// The Welch factors describe measured bin second moments, not deterministic
/// waveform phases across bins. Overlap supplies no independence, confidence
/// interval or stationarity guarantee. Finite-record leakage, source aliasing,
/// spatial quadrature and the propagation owner's physical scope remain.
pub struct RecordedQuadrupoleSpectrum {
    source_positions_m: Vec<[f64; 3]>,
    clock: QuadrupoleSampleClock,
    spectrum: WelchCrossSpectrum,
}

impl RecordedQuadrupoleSpectrum {
    pub const fn clock(&self) -> QuadrupoleSampleClock {
        self.clock
    }

    pub fn source_positions_m(&self) -> &[[f64; 3]] {
        &self.source_positions_m
    }

    pub fn spectrum(&self) -> &WelchCrossSpectrum {
        &self.spectrum
    }

    /// Prepare the original pressure/vector-velocity transfer at an interior
    /// recorded bin. Source positions cannot be silently reassigned.
    pub fn prepare_observer(
        &self,
        acoustic: FluidAcousticState,
        observer_m: [f64; 3],
        bin: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedQuadrupoleObserver, AeroacError> {
        poll(&mut cancelled)?;
        let frequency_hz = self
            .spectrum
            .frequency_hz(bin)
            .map_err(AeroacError::Welch)?;
        self.spectrum
            .peak_bin_hz(bin, frequency_hz)
            .map_err(AeroacError::Welch)?;
        PreparedQuadrupoleObserver::new_cancellable(
            acoustic,
            frequency_hz,
            observer_m,
            &self.source_positions_m,
            self.source_positions_m.len(),
            cancelled,
        )
    }

    /// Bind the same recorded source layout and bin to the original molecular
    /// gas. Absorption and dispersion affect every tensor field coefficient;
    /// the retained Welch factors, source units and one-sided scale are unchanged.
    pub fn prepare_molecular_observer(
        &self,
        source: RelaxationState,
        observer_m: [f64; 3],
        bin: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedQuadrupoleObserver, MaterialAeroacError> {
        poll(&mut cancelled)?;
        let frequency_hz = self
            .spectrum
            .frequency_hz(bin)
            .map_err(AeroacError::Welch)?;
        self.spectrum
            .peak_bin_hz(bin, frequency_hz)
            .map_err(AeroacError::Welch)?;
        PreparedQuadrupoleObserver::new_molecular_cancellable(
            source,
            frequency_hz,
            observer_m,
            &self.source_positions_m,
            self.source_positions_m.len(),
            cancelled,
        )
    }

    /// Coherently propagate each recorded factor through its bound geometry,
    /// then form bin-integrated pressure/velocity moments and signed intensity.
    pub fn evaluate_bin(
        &self,
        observer: &PreparedQuadrupoleObserver,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        poll(&mut cancelled)?;
        if observer.source_count() != self.source_positions_m.len() {
            return Err(invalid("quadrupole observer has a different source layout"));
        }
        for (i, &position) in self.source_positions_m.iter().enumerate() {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            if observer.source_position_m(i) != Some(position) {
                return Err(invalid(
                    "quadrupole observer changed recorded source geometry",
                ));
            }
        }
        observer.evaluate_welch_bin(&self.spectrum, bin, budget, cancelled)
    }
}

fn invalid(what: &'static str) -> AeroacError {
    AeroacError::InvalidParameter { what }
}

#[cfg(test)]
mod tests;
