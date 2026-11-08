//! One synchronized record for disjoint compact and volume-stress sources.
//!
//! Physical `[Q,Fx,Fy,Fz]` blocks precede `[xx,yy,zz,xy,xz,yz]` stress blocks.
//! One original Welch owner supplies every window, segment and normalization.
//! Its factors enter the existing mixed pressure/vector-velocity transfer
//! before moments are formed, retaining cross-source and cross-mechanism work.
//! Manual channels remain caller-owned: the caller establishes a common
//! physical history and disjoint source mechanisms. In particular, enclosed
//! volume stresses must not be added to a permeable FW-H surface a second time.

use super::SourceSampleClock;
use crate::AeroacError;
use crate::curle3d::array::{PreparedCompactSourceObserver, plane};
use crate::curle3d::quadrupole::PreparedQuadrupoleObserver;
use crate::curle3d::statistics::invalid;
use crate::curle3d::{CoherenceBudget, CompactFieldStatistics};
use crate::material::MaterialAeroacError;
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchStream};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::relaxation::RelaxationState;
use fs_math::c64::C64;

/// Bounds both original source layouts and the SINGLE joint transform.
#[derive(Clone, Copy, Debug)]
pub struct JointSourceRecordingBudget {
    pub maximum_compact_sources: usize,
    pub maximum_quadrupole_sources: usize,
    /// Includes both retained position arrays and the original Welch buffers.
    /// Caller histories and allocator metadata are excluded. Transform work
    /// counts all 4*S + 6*V channels, with no separate factor normalization.
    pub welch: WelchBudget,
}

/// Fixed geometry, reference and clock for two coherent physical mechanisms.
///
/// Each group has at least one source; an absent mechanism is an explicit zero
/// sample. For a compact-only record use the existing CompactSourceWelchStream.
/// Pushes allocate nothing and are atomic even if cancelled during an FFT.
pub struct JointSourceWelchStream {
    acoustic: FluidAcousticState,
    compact_positions_m: Vec<[f64; 3]>,
    quadrupole_positions_m: Vec<[f64; 3]>,
    clock: SourceSampleClock,
    stream: WelchStream,
}

impl JointSourceWelchStream {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        acoustic: FluidAcousticState,
        compact_positions_m: &[[f64; 3]],
        quadrupole_positions_m: &[[f64; 3]],
        clock: SourceSampleClock,
        expected_samples: usize,
        config: WelchConfig,
        budget: JointSourceRecordingBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        plane::poll(&mut cancelled)?;
        let (compact, stress) = (compact_positions_m.len(), quadrupole_positions_m.len());
        if compact == 0
            || stress == 0
            || compact > budget.maximum_compact_sources
            || stress > budget.maximum_quadrupole_sources
            || expected_samples == 0
        {
            return Err(invalid(
                "joint recording needs bounded source groups and samples",
            ));
        }
        clock.validate_last_index(expected_samples - 1)?;
        let channels = compact
            .checked_mul(4)
            .and_then(|n| stress.checked_mul(6).and_then(|m| n.checked_add(m)))
            .ok_or_else(|| invalid("joint recording channel count overflow"))?;
        let geometry_bytes = compact
            .checked_add(stress)
            .and_then(|n| n.checked_mul(size_of::<[f64; 3]>()))
            .ok_or_else(|| invalid("joint recording geometry extent overflow"))?;
        let mut transform_budget = budget.welch;
        transform_budget.maximum_workspace_bytes = transform_budget
            .maximum_workspace_bytes
            .checked_sub(geometry_bytes)
            .ok_or_else(|| invalid("joint geometry exceeds recording memory budget"))?;
        let compact_positions_m = copy_positions(compact_positions_m, &mut cancelled)?;
        let quadrupole_positions_m = copy_positions(quadrupole_positions_m, &mut cancelled)?;
        let stream = WelchStream::new(
            channels,
            expected_samples,
            clock.timestep_s(),
            config,
            transform_budget,
            &mut cancelled,
        )
        .map_err(AeroacError::Welch)?;
        plane::poll(&mut cancelled)?;
        Ok(Self {
            acoustic,
            compact_positions_m,
            quadrupole_positions_m,
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
    pub fn compact_source_positions_m(&self) -> &[[f64; 3]] {
        &self.compact_positions_m
    }
    pub fn quadrupole_source_positions_m(&self) -> &[[f64; 3]] {
        &self.quadrupole_positions_m
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

    /// Append one simultaneous real physical sample. Q is m³/s, force ON the
    /// fluid is N, and symmetric stress moments are N m with physical (not
    /// pre-doubled) shear entries. No area, volume, density, derivative or
    /// source phase is added. The source order and geometry remain immutable.
    /// Any refusal leaves both groups' accepted clock and factors unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn push_sample(
        &mut self,
        compact_positions_m: &[[f64; 3]],
        quadrupole_positions_m: &[[f64; 3]],
        sample_time_s: f64,
        volume_force: &[[f64; 4]],
        stress_moments_n_m: &[[f64; 6]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        plane::poll(&mut cancelled)?;
        if volume_force.len() != self.compact_positions_m.len()
            || stress_moments_n_m.len() != self.quadrupole_positions_m.len()
            || self.remaining_sample_count() == 0
        {
            return Err(invalid("joint recording sample shape or declared length"));
        }
        self.clock
            .admit(self.accepted_sample_count(), sample_time_s)?;
        same_positions(
            compact_positions_m,
            &self.compact_positions_m,
            &mut cancelled,
        )?;
        same_positions(
            quadrupole_positions_m,
            &self.quadrupole_positions_m,
            &mut cancelled,
        )?;
        let compact_channels = 4 * self.compact_positions_m.len();
        self.stream
            .push_with(
                compact_channels + 6 * self.quadrupole_positions_m.len(),
                1,
                |channel, _| {
                    if channel < compact_channels {
                        volume_force[channel / 4][channel % 4]
                    } else {
                        let stress = channel - compact_channels;
                        stress_moments_n_m[stress / 6][stress % 6]
                    }
                },
                &mut cancelled,
            )
            .map_err(AeroacError::Welch)
    }

    /// Publish only a complete declared record, without copying its factors.
    pub fn finish(self) -> Result<RecordedJointSourceSpectrum, AeroacError> {
        Ok(RecordedJointSourceSpectrum {
            acoustic: self.acoustic,
            compact_positions_m: self.compact_positions_m,
            quadrupole_positions_m: self.quadrupole_positions_m,
            clock: self.clock,
            spectrum: self.stream.finish().map_err(AeroacError::Welch)?,
        })
    }
}

/// Original mixed source factors and their two immutable physical layouts.
///
/// This is a bin-integrated second-moment observation, not a reconstruction of
/// cross-frequency phases. Finite-record leakage, source/volume resolution,
/// aliasing and the native propagation domain remain. Overlap implies neither
/// independent observations nor a confidence interval.
pub struct RecordedJointSourceSpectrum {
    acoustic: FluidAcousticState,
    compact_positions_m: Vec<[f64; 3]>,
    quadrupole_positions_m: Vec<[f64; 3]>,
    clock: SourceSampleClock,
    spectrum: WelchCrossSpectrum,
}

impl RecordedJointSourceSpectrum {
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }
    pub const fn clock(&self) -> SourceSampleClock {
        self.clock
    }
    pub fn compact_source_positions_m(&self) -> &[[f64; 3]] {
        &self.compact_positions_m
    }
    pub fn quadrupole_source_positions_m(&self) -> &[[f64; 3]] {
        &self.quadrupole_positions_m
    }
    pub fn spectrum(&self) -> &WelchCrossSpectrum {
        &self.spectrum
    }

    /// Prepare both native full-vector transfers at the same interior bin.
    /// Callers may add matching ideal planes/wedges to both owners; mixed
    /// evaluation refuses differing boundary domains or material laws.
    pub fn prepare_observers(
        &self,
        observer_m: [f64; 3],
        bin: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(PreparedCompactSourceObserver, PreparedQuadrupoleObserver), AeroacError> {
        plane::poll(&mut cancelled)?;
        let frequency = self.bin_frequency(bin)?;
        let compact = PreparedCompactSourceObserver::prepare_cancellable(
            self.acoustic,
            frequency,
            observer_m,
            &self.compact_positions_m,
            self.compact_positions_m.len(),
            C64::from_re(core::f64::consts::TAU * frequency / self.acoustic.sound_speed_m_s()),
            None,
            &mut cancelled,
        )?;
        let stress = PreparedQuadrupoleObserver::new_cancellable(
            self.acoustic,
            frequency,
            observer_m,
            &self.quadrupole_positions_m,
            self.quadrupole_positions_m.len(),
            &mut cancelled,
        )?;
        Ok((compact, stress))
    }

    /// Both mechanisms retain the SAME original molecular bulk response.
    /// Equilibrium rho/c must match the record; no independent attenuation,
    /// extra phase or plane-wave velocity is authored by this adapter.
    pub fn prepare_molecular_observers(
        &self,
        source: RelaxationState,
        observer_m: [f64; 3],
        bin: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(PreparedCompactSourceObserver, PreparedQuadrupoleObserver), MaterialAeroacError>
    {
        plane::poll(&mut cancelled)?;
        if source
            .gas()
            .acoustic_state()
            .map_err(MaterialAeroacError::Material)?
            != self.acoustic
        {
            return Err(invalid("joint record and molecular acoustic references differ").into());
        }
        let frequency = self.bin_frequency(bin)?;
        let k = source
            .relaxing_wavenumber(core::f64::consts::TAU * frequency)
            .map_err(MaterialAeroacError::Material)?;
        let compact = PreparedCompactSourceObserver::prepare_cancellable(
            self.acoustic,
            frequency,
            observer_m,
            &self.compact_positions_m,
            self.compact_positions_m.len(),
            k,
            Some(source.clone()),
            &mut cancelled,
        )?;
        let stress = PreparedQuadrupoleObserver::new_molecular_cancellable(
            source,
            frequency,
            observer_m,
            &self.quadrupole_positions_m,
            self.quadrupole_positions_m.len(),
            &mut cancelled,
        )?;
        Ok((compact, stress))
    }

    /// Propagate each original factor coherently through BOTH mechanisms,
    /// then form pressure/velocity covariance and signed active/reactive flux.
    /// The native mixed operator applies sqrt(2 df) exactly once and admits
    /// the combined source/factor work, frequency and propagation domain.
    pub fn evaluate_bin(
        &self,
        compact: &PreparedCompactSourceObserver,
        stress: &PreparedQuadrupoleObserver,
        bin: usize,
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        plane::poll(&mut cancelled)?;
        if compact.acoustic() != self.acoustic
            || stress.acoustic() != self.acoustic
            || compact.source_count() != self.compact_positions_m.len()
            || stress.source_count() != self.quadrupole_positions_m.len()
        {
            return Err(invalid(
                "joint observer changed source counts or acoustic reference",
            ));
        }
        for (index, &position) in self.compact_positions_m.iter().enumerate() {
            if index.is_multiple_of(256) {
                plane::poll(&mut cancelled)?;
            }
            if compact.source_position_m(index) != Some(position) {
                return Err(invalid(
                    "joint compact observer changed source geometry or order",
                ));
            }
        }
        for (index, &position) in self.quadrupole_positions_m.iter().enumerate() {
            if index.is_multiple_of(256) {
                plane::poll(&mut cancelled)?;
            }
            if stress.source_position_m(index) != Some(position) {
                return Err(invalid(
                    "joint stress observer changed source geometry or order",
                ));
            }
        }
        compact.evaluate_welch_bin_with_quadrupoles(stress, &self.spectrum, bin, budget, cancelled)
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

fn copy_positions(
    positions: &[[f64; 3]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<[f64; 3]>, AeroacError> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(positions.len())
        .map_err(|_| invalid("joint source geometry allocation"))?;
    for (index, &position) in positions.iter().enumerate() {
        if index.is_multiple_of(256) {
            plane::poll(cancelled)?;
        }
        if position.iter().any(|x| !x.is_finite()) {
            return Err(invalid("joint source geometry is nonfinite"));
        }
        output.push(position);
    }
    Ok(output)
}

fn same_positions(
    actual: &[[f64; 3]],
    expected: &[[f64; 3]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AeroacError> {
    if actual.len() != expected.len() {
        return Err(invalid("joint recording source geometry count changed"));
    }
    for (index, (a, b)) in actual.iter().zip(expected).enumerate() {
        if index.is_multiple_of(256) {
            plane::poll(cancelled)?;
        }
        if a != b {
            return Err(invalid("joint recording source geometry or order changed"));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
