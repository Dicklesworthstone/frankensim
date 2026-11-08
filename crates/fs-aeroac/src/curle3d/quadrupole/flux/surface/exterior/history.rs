//! Native FW-H surface plus exterior-volume HISTORY from one accepted flow.
//!
//! The original extractors read the same immutable Euler/mixture state. Their
//! physical Q/F/T samples, clocks, frame and source ordering remain unchanged.
//! Preparation explicitly chooses an offline clamped C2 reconstruction; causal
//! observation declares zero past and uses the original molecular Q/F/T owner.
//! This does not infer derivatives from an unqualified noisy signal, reconstruct
//! Welch phases, resample CFD, advance flow or synthesize an unrecorded tail.

use super::{
    CartesianGrid, CompressibleVolumeConfig3, CompressibleVolumeError3, Euler3d,
    ExteriorSurfaceRecordingError, MixtureEuler3d, PermeableBox3, PermeableSurfaceControls3,
    PreparedCompressibleVolume3, PreparedPermeableSurface3, SourceSampleClock, admit_exterior,
    admit_exterior_receiver, invalid, poll,
};
use crate::AeroacError;
use crate::curle3d::quadrupole::transient::history::QuadrupoleHistoryEnds;
use crate::curle3d::quadrupole::transient::molecular::sources::history::{
    AcousticSourceHistoryBudget, AcousticSourceHistoryRecorder, CompactHistoryEnds,
    PreparedAcousticSourceHistory,
};
use crate::curle3d::quadrupole::transient::molecular::sources::{
    MolecularSourceRecord, MolecularTransientBudget, MolecularTransientConfig,
    PreparedMolecularSourceObserver,
};
use crate::material::MaterialAeroacError;
use fs_material::relaxation::RelaxationState;

/// Peak native storage is bounded by the SUM of these independent caps.
/// The recording cap includes simultaneous raw and prepared histories. Native
/// material/flow metadata and the separately admitted propagation are excluded.
#[derive(Clone, Copy, Debug)]
pub struct ExteriorSurfaceHistoryBudget {
    pub surface: PermeableSurfaceControls3,
    pub maximum_volume_workspace_bytes: usize,
    pub recording: AcousticSourceHistoryBudget,
}

/// Both source samplers remain private; independently sampled same-time caches
/// cannot substitute for the ONE accepted flow borrow used by each sample.
/// Constructors record the current accepted state exactly once. Later samples
/// must satisfy the same original uniform clock without hidden resampling.
pub struct ExteriorSurfaceHistoryRecorder {
    surface: PreparedPermeableSurface3,
    volume: PreparedCompressibleVolume3,
    recording: AcousticSourceHistoryRecorder,
}

impl ExteriorSurfaceHistoryRecorder {
    #[allow(clippy::too_many_arguments)]
    pub fn from_frozen(
        flow: &Euler3d,
        control_box: PermeableBox3,
        exterior_cells: &[usize],
        source_config: CompressibleVolumeConfig3,
        clock: SourceSampleClock,
        expected_samples: usize,
        budget: ExteriorSurfaceHistoryBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ExteriorSurfaceRecordingError> {
        poll(&mut cancelled)?;
        clock.admit(0, flow.time().value())?;
        admit_exterior(
            flow.grid(),
            control_box,
            exterior_cells,
            source_config,
            budget.recording.maximum_compact_sources,
            budget.recording.maximum_quadrupole_sources,
            budget.maximum_volume_workspace_bytes,
            &mut cancelled,
        )?;
        let surface = PreparedPermeableSurface3::from_frozen(
            flow,
            control_box,
            source_config,
            budget.surface,
            &mut cancelled,
        )?;
        let volume = PreparedCompressibleVolume3::from_frozen(
            flow,
            exterior_cells,
            source_config,
            budget.maximum_volume_workspace_bytes,
            &mut cancelled,
        )?;
        Self::from_sources(
            surface,
            volume,
            clock,
            expected_samples,
            budget.recording,
            cancelled,
        )
    }

    /// The complete original species/caloric owner and explicit ambient
    /// composition are retained by the native surface. Current cells retain
    /// their own compositions in BOTH source extractions from the same flow.
    #[allow(clippy::too_many_arguments)]
    pub fn from_mixture(
        flow: &MixtureEuler3d,
        reference_partial_densities_kg_m3: &[f64],
        control_box: PermeableBox3,
        exterior_cells: &[usize],
        source_config: CompressibleVolumeConfig3,
        clock: SourceSampleClock,
        expected_samples: usize,
        budget: ExteriorSurfaceHistoryBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ExteriorSurfaceRecordingError> {
        poll(&mut cancelled)?;
        clock.admit(0, flow.time().value())?;
        admit_exterior(
            flow.grid(),
            control_box,
            exterior_cells,
            source_config,
            budget.recording.maximum_compact_sources,
            budget.recording.maximum_quadrupole_sources,
            budget.maximum_volume_workspace_bytes,
            &mut cancelled,
        )?;
        let surface = PreparedPermeableSurface3::from_mixture(
            flow,
            reference_partial_densities_kg_m3,
            control_box,
            source_config,
            budget.surface,
            &mut cancelled,
        )?;
        let volume = PreparedCompressibleVolume3::from_mixture(
            flow,
            exterior_cells,
            source_config,
            budget.maximum_volume_workspace_bytes,
            &mut cancelled,
        )?;
        Self::from_sources(
            surface,
            volume,
            clock,
            expected_samples,
            budget.recording,
            cancelled,
        )
    }

    fn from_sources(
        surface: PreparedPermeableSurface3,
        volume: PreparedCompressibleVolume3,
        clock: SourceSampleClock,
        expected_samples: usize,
        budget: AcousticSourceHistoryBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ExteriorSurfaceRecordingError> {
        let recording = AcousticSourceHistoryRecorder::new(
            surface.config.acoustic,
            &surface.positions,
            &volume.positions,
            clock,
            expected_samples,
            budget,
            &mut cancelled,
        )?;
        let mut result = Self {
            surface,
            volume,
            recording,
        };
        result.push_accepted(&mut cancelled)?;
        // Nothing fallible follows the single successful recorder commit.
        Ok(result)
    }

    /// Read one accepted native state for both source mechanisms. No advance,
    /// resampling or warm allocation; every refusal leaves both source caches
    /// and the joint recording at their preceding accepted sample.
    pub fn sample_frozen(
        &mut self,
        flow: &Euler3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ExteriorSurfaceRecordingError> {
        self.admit_time(flow.time().value(), &mut cancelled)?;
        let previous = (self.surface.time_s, self.volume.time_s);
        self.surface.sample_frozen(flow, &mut cancelled)?;
        let volume_result = self.volume.sample_frozen(flow, &mut cancelled);
        self.finish_sample(previous, volume_result, &mut cancelled)
    }

    pub fn sample_mixture(
        &mut self,
        flow: &MixtureEuler3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ExteriorSurfaceRecordingError> {
        self.admit_time(flow.time().value(), &mut cancelled)?;
        let previous = (self.surface.time_s, self.volume.time_s);
        self.surface.sample_mixture(flow, &mut cancelled)?;
        let volume_result = self.volume.sample_mixture(flow, &mut cancelled);
        self.finish_sample(previous, volume_result, &mut cancelled)
    }

    fn admit_time(
        &self,
        time_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ExteriorSurfaceRecordingError> {
        poll(cancelled)?;
        if self.remaining_sample_count() == 0 {
            return Err(invalid("joint FW-H declared record is already complete").into());
        }
        self.clock().admit(self.accepted_sample_count(), time_s)?;
        Ok(())
    }

    // The surface succeeded before entry. A failed native sampler has NOT
    // committed: its dirty trial must never be swapped into accepted storage.
    fn finish_sample(
        &mut self,
        previous: (Option<f64>, Option<f64>),
        volume_result: Result<(), CompressibleVolumeError3>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ExteriorSurfaceRecordingError> {
        if let Err(error) = volume_result {
            self.restore_surface(previous.0);
            return Err(error.into());
        }
        if let Err(error) = self.push_accepted(cancelled) {
            self.restore_surface(previous.0);
            core::mem::swap(&mut self.volume.moments, &mut self.volume.trial);
            self.volume.time_s = previous.1;
            return Err(error);
        }
        // After the raw sample succeeds, no further poll or refusal can invalidate it.
        Ok(())
    }

    fn restore_surface(&mut self, previous: Option<f64>) {
        core::mem::swap(&mut self.surface.sources, &mut self.surface.trial);
        self.surface.time_s = previous;
    }

    fn push_accepted(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ExteriorSurfaceRecordingError> {
        if self.surface.time_s != self.volume.time_s {
            return Err(invalid("joint FW-H native sample times differ").into());
        }
        self.recording.push_sample(
            &self.surface.positions,
            &self.volume.positions,
            self.surface.sample_time_s(),
            &self.surface.sources,
            &self.volume.moments,
            cancelled,
        )?;
        Ok(())
    }

    pub fn surface(&self) -> &PreparedPermeableSurface3 {
        &self.surface
    }
    pub fn volume(&self) -> &PreparedCompressibleVolume3 {
        &self.volume
    }
    pub const fn clock(&self) -> SourceSampleClock {
        self.recording.clock()
    }
    pub fn accepted_sample_count(&self) -> usize {
        self.recording.accepted_sample_count()
    }
    pub fn remaining_sample_count(&self) -> usize {
        self.recording.remaining_sample_count()
    }
    /// Read-only accepted raw samples; this does not expose a mutable sampler.
    pub fn recording(&self) -> &AcousticSourceHistoryRecorder {
        &self.recording
    }

    /// Prepare one complete original history without consuming the recording.
    /// Endpoint physics and the below-Nyquist bandwidth declaration are explicit;
    /// a cancelled or refused reconstruction can be retried without resampling.
    pub fn prepare_clamped_cubic(
        &self,
        compact_ends: &[CompactHistoryEnds],
        stress_ends: &[QuadrupoleHistoryEnds],
        maximum_source_frequency_hz: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RecordedExteriorSurfaceHistory, ExteriorSurfaceRecordingError> {
        let history = self.recording.prepare_clamped_cubic(
            compact_ends,
            stress_ends,
            maximum_source_frequency_hz,
            cancelled,
        )?;
        Ok(RecordedExteriorSurfaceHistory {
            grid: self.surface.grid,
            control_box: self.surface.control_box,
            source_config: self.surface.config,
            history,
        })
    }
}

/// Complete transient source record retaining the FW-H exterior receiver domain.
/// Enclosed volume cells and receivers on/inside the closed surface are refused.
/// Selected exterior extent is finite; inviscid, stationary-source and source
/// reconstruction limits remain. Molecular output is an initially-resting finite
/// prefix, not a claim of complete later echoes or a calibrated audio render.
pub struct RecordedExteriorSurfaceHistory {
    grid: CartesianGrid,
    control_box: PermeableBox3,
    source_config: CompressibleVolumeConfig3,
    history: PreparedAcousticSourceHistory,
}

impl RecordedExteriorSurfaceHistory {
    pub const fn clock(&self) -> SourceSampleClock {
        self.history.clock()
    }
    pub const fn control_box(&self) -> PermeableBox3 {
        self.control_box
    }
    pub const fn source_config(&self) -> CompressibleVolumeConfig3 {
        self.source_config
    }
    pub const fn sample_count(&self) -> usize {
        self.history.sample_count()
    }
    pub fn source_positions_m(&self) -> &[[f64; 3]] {
        self.history.source_positions_m()
    }
    pub fn compact_source_count(&self) -> usize {
        self.history.compact_source_count()
    }
    pub fn quadrupole_source_count(&self) -> usize {
        self.history.quadrupole_source_count()
    }
    pub const fn declared_maximum_source_frequency_hz(&self) -> f64 {
        self.history.declared_maximum_source_frequency_hz()
    }
    pub const fn workspace_bytes(&self) -> usize {
        self.history.workspace_bytes()
    }

    pub fn prepare_molecular_observer(
        &self,
        molecular: RelaxationState,
        receiver_m: [f64; 3],
        config: MolecularTransientConfig,
        budget: MolecularTransientBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularSourceObserver, MaterialAeroacError> {
        admit_exterior_receiver(&self.grid, self.control_box, self.source_config, receiver_m)?;
        self.history
            .prepare_molecular_observer(molecular, receiver_m, config, budget, cancelled)
    }

    /// Enforce the original FW-H domain, source identity, clock and material
    /// reference before one coherent complete pressure/vector-velocity record.
    /// Calling explicitly declares zero prehistory and unexcited material; the
    /// original propagation owner refuses incompatible onset or an uncovered
    /// output prefix. No source attenuation or propagation delay is duplicated.
    pub fn observe_initially_resting(
        &self,
        observer: &PreparedMolecularSourceObserver,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularSourceRecord, AeroacError> {
        admit_exterior_receiver(
            &self.grid,
            self.control_box,
            self.source_config,
            observer.paths().observer_position_m(),
        )?;
        self.history.observe_initially_resting(observer, cancelled)
    }
}

#[cfg(test)]
mod tests;
