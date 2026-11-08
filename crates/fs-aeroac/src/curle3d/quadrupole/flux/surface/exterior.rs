//! Coherent FW-H surface and selected EXTERIOR Lighthill volume acquisition.
//!
//! Both native extractors read ONE borrowed accepted compressible-flow state.
//! Their physical mass/momentum and stress moments enter one joint Welch push;
//! no source kernel, derivative, density, area or cell-volume factor is added.
//! The native samplers remain private. A rejected observation restores their
//! accepted buffers/times and the recorder can retry the same physical sample.
//!
//! Only original grid cells outside the control box are admitted. An exterior
//! cell sharing a surface face remains valid: its participation in centered
//! face reconstruction does not make its volume an enclosed source. Selected
//! exterior cells do not complete an infinite volume; extent and spatial
//! refinement remain necessary. No viscosity, SGS stress, body-force dipoles,
//! moving geometry, CFD feedback or automatic sound-validity claim is added.

use super::super::PreparedCompressibleVolume3;
use super::{
    CartesianGrid, CompressibleVolumeConfig3, CompressibleVolumeError3, Euler3d, MixtureEuler3d,
    PermeableBox3, PermeableSurfaceControls3, PreparedPermeableSurface3, invalid, poll,
};
use crate::AeroacError;
use crate::curle3d::array::recorded::SourceSampleClock;
use crate::curle3d::array::recorded::joint::{
    JointRadiatedPowerBand, JointRadiatedPowerBandBudget, JointSourceRecordingBudget,
    JointSourceWelchStream, PreparedJointRadiatedPower, RecordedJointSourceSpectrum,
};
use crate::curle3d::power::RadiatedPowerBudget;
use crate::curle3d::quadrupole::PreparedQuadrupoleObserver;
use crate::curle3d::{CoherenceBudget, CompactFieldStatistics, PreparedCompactSourceObserver};
use crate::material::MaterialAeroacError;
use fs_fft::welch::{WelchConfig, WelchCrossSpectrum};
use fs_material::relaxation::RelaxationState;
use fs_math::dd::Dd;

/// Accepted real source histories for complete causal molecular observation.
pub mod history;

/// Independent, explicit caps for both samplers and their joint spectrum.
/// Peak numerical storage is bounded by the SUM of these three byte caps;
/// gas/species metadata, caller flow storage and allocator overhead are excluded.
#[derive(Clone, Copy, Debug)]
pub struct ExteriorSurfaceRecordingBudget {
    pub surface: PermeableSurfaceControls3,
    pub maximum_volume_workspace_bytes: usize,
    /// Source counts also bound extraction work before allocating either cache.
    /// One Welch transform budget counts every surface and volume channel.
    pub recording: JointSourceRecordingBudget,
}

/// Preserve original flow/EOS errors separately from acoustic observation errors.
#[derive(Debug)]
pub enum ExteriorSurfaceRecordingError {
    Source(CompressibleVolumeError3),
    Acoustic(AeroacError),
}
impl core::fmt::Display for ExteriorSurfaceRecordingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(error) => write!(f, "joint FW-H source: {error}"),
            Self::Acoustic(error) => write!(f, "joint FW-H recording: {error}"),
        }
    }
}
impl std::error::Error for ExteriorSurfaceRecordingError {}
impl From<CompressibleVolumeError3> for ExteriorSurfaceRecordingError {
    fn from(error: CompressibleVolumeError3) -> Self {
        Self::Source(error)
    }
}
impl From<AeroacError> for ExteriorSurfaceRecordingError {
    fn from(error: AeroacError) -> Self {
        Self::Acoustic(error)
    }
}

/// Native stationary surface/exterior-volume record from one physical solver.
///
/// Constructors include the current accepted flow sample exactly once. Later
/// calls observe the next declared uniform-clock sample without advancing CFD.
/// Full original gas, grid, boundary topology, frame/reference and source order
/// stay bound through the native owners; independently authored same-time
/// source snapshots cannot be substituted into this primary path.
pub struct ExteriorSurfaceWelchStream {
    surface: PreparedPermeableSurface3,
    volume: PreparedCompressibleVolume3,
    recording: JointSourceWelchStream,
}

impl ExteriorSurfaceWelchStream {
    #[allow(clippy::too_many_arguments)]
    pub fn from_frozen(
        flow: &Euler3d,
        control_box: PermeableBox3,
        exterior_cells: &[usize],
        source_config: CompressibleVolumeConfig3,
        clock: SourceSampleClock,
        expected_samples: usize,
        config: WelchConfig,
        budget: ExteriorSurfaceRecordingBudget,
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
            config,
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
        config: WelchConfig,
        budget: ExteriorSurfaceRecordingBudget,
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
            config,
            budget.recording,
            cancelled,
        )
    }

    fn from_sources(
        surface: PreparedPermeableSurface3,
        volume: PreparedCompressibleVolume3,
        clock: SourceSampleClock,
        expected_samples: usize,
        config: WelchConfig,
        budget: JointSourceRecordingBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ExteriorSurfaceRecordingError> {
        let recording = JointSourceWelchStream::new(
            surface.config.acoustic,
            &surface.positions,
            &volume.positions,
            clock,
            expected_samples,
            config,
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
        // After Welch succeeds, no further poll or refusal can invalidate it.
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
    pub fn completed_segment_count(&self) -> usize {
        self.recording.completed_segment_count()
    }

    /// Retain the FW-H observation domain beside the complete joint factors.
    /// Check remaining_sample_count before consuming an incomplete record.
    pub fn finish(self) -> Result<RecordedExteriorSurfaceSpectrum, ExteriorSurfaceRecordingError> {
        Ok(RecordedExteriorSurfaceSpectrum {
            grid: self.surface.grid,
            control_box: self.surface.control_box,
            source_config: self.surface.config,
            recorded: self.recording.finish()?,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn admit_exterior(
    grid: &CartesianGrid,
    control_box: PermeableBox3,
    cells: &[usize],
    config: CompressibleVolumeConfig3,
    maximum_compact_sources: usize,
    maximum_quadrupole_sources: usize,
    maximum_volume_workspace_bytes: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), CompressibleVolumeError3> {
    poll(cancelled)?;
    if control_box.face_count(grid)? > maximum_compact_sources
        || cells.len() > maximum_quadrupole_sources
    {
        return Err(invalid("joint FW-H source counts exceed recording budget"));
    }
    // Reuse native bounds, increasing-order and physical-frame admission.
    PreparedCompressibleVolume3::admit(
        grid,
        cells,
        config,
        maximum_volume_workspace_bytes,
        cancelled,
    )?;
    let [nx, ny, _] = grid.shape();
    for (index, &cell) in cells.iter().enumerate() {
        if index.is_multiple_of(64) {
            poll(cancelled)?;
        }
        let xyz = [cell % nx, (cell / nx) % ny, cell / nx / ny];
        if (0..3).all(|a| control_box.lower[a] <= xyz[a] && xyz[a] < control_box.upper[a]) {
            return Err(invalid(
                "joint FW-H volume contains an enclosed control-box cell",
            ));
        }
    }
    Ok(())
}

/// Complete acquired joint factors with the original FW-H exterior domain.
///
/// Pressure interpretation requires a receiver OUTSIDE the closed control box.
/// These preparation/evaluation methods enforce that domain; raw spectrum
/// access is for physical-channel analysis, not an interior-pressure claim.
/// Finite selected exterior extent and all inviscid source limits remain.
pub struct RecordedExteriorSurfaceSpectrum {
    grid: CartesianGrid,
    control_box: PermeableBox3,
    source_config: CompressibleVolumeConfig3,
    recorded: RecordedJointSourceSpectrum,
}

impl RecordedExteriorSurfaceSpectrum {
    pub const fn clock(&self) -> SourceSampleClock {
        self.recorded.clock()
    }
    pub fn spectrum(&self) -> &WelchCrossSpectrum {
        self.recorded.spectrum()
    }
    pub const fn control_box(&self) -> PermeableBox3 {
        self.control_box
    }
    pub const fn source_config(&self) -> CompressibleVolumeConfig3 {
        self.source_config
    }
    pub fn compact_source_positions_m(&self) -> &[[f64; 3]] {
        self.recorded.compact_source_positions_m()
    }
    pub fn quadrupole_source_positions_m(&self) -> &[[f64; 3]] {
        self.recorded.quadrupole_source_positions_m()
    }

    /// Full-sphere lossless power from this SAME native surface/volume record.
    /// The original source partition, positions, clock and measured coherence
    /// enter one analytic multipole power operator. No microphone or separate
    /// source reconstruction is needed; the finite exterior-volume scope remains.
    pub fn prepare_radiated_power(
        &self,
        bin: usize,
        budget: RadiatedPowerBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedJointRadiatedPower<'_>, AeroacError> {
        self.recorded.prepare_radiated_power(bin, budget, cancelled)
    }

    /// Integrate complete native interior bins with every surface/stress cross term.
    /// The joint owner bounds the complete sweep and returns W/Hz rows plus W.
    /// This is a stationary spectral observation, without absorption or inferred
    /// contributions from unselected exterior cells or unresolved source physics.
    pub fn radiated_power_band(
        &self,
        bins: core::ops::Range<usize>,
        budget: JointRadiatedPowerBandBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<JointRadiatedPowerBand, AeroacError> {
        self.recorded.radiated_power_band(bins, budget, cancelled)
    }

    pub fn prepare_observers(
        &self,
        observer_m: [f64; 3],
        bin: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(PreparedCompactSourceObserver, PreparedQuadrupoleObserver), AeroacError> {
        self.admit_receiver(observer_m)?;
        self.recorded.prepare_observers(observer_m, bin, cancelled)
    }

    pub fn prepare_molecular_observers(
        &self,
        source: RelaxationState,
        observer_m: [f64; 3],
        bin: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(PreparedCompactSourceObserver, PreparedQuadrupoleObserver), MaterialAeroacError>
    {
        self.admit_receiver(observer_m)?;
        self.recorded
            .prepare_molecular_observers(source, observer_m, bin, cancelled)
    }

    pub fn evaluate_bin(
        &self,
        compact: &PreparedCompactSourceObserver,
        stress: &PreparedQuadrupoleObserver,
        bin: usize,
        budget: CoherenceBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.admit_receiver(compact.observer_position_m())?;
        // The native mixed owner also requires exactly equal receiver positions.
        self.recorded
            .evaluate_bin(compact, stress, bin, budget, cancelled)
    }

    fn admit_receiver(&self, point: [f64; 3]) -> Result<(), AeroacError> {
        admit_exterior_receiver(&self.grid, self.control_box, self.source_config, point)
    }
}

// The same closed-box receiver domain applies to spectral and transient paths.
fn admit_exterior_receiver(
    grid: &CartesianGrid,
    control_box: PermeableBox3,
    source_config: CompressibleVolumeConfig3,
    point: [f64; 3],
) -> Result<(), AeroacError> {
    let invalid = || AeroacError::InvalidParameter {
        what: "FW-H pressure receiver must lie outside the closed control box",
    };
    if point.iter().any(|x| !x.is_finite()) {
        return Err(invalid());
    }
    let delta: [Dd; 3] =
        core::array::from_fn(|a| Dd::from_f64(point[a]) - Dd::from_f64(source_config.origin_m[a]));
    let columns: [[Dd; 3]; 3] = core::array::from_fn(|a| {
        core::array::from_fn(|b| Dd::from_f64(source_config.grid_to_acoustic[b][a]))
    });
    let spacing = grid.spacing();
    let scale = point
        .iter()
        .chain(&source_config.origin_m)
        .map(|x| x.abs())
        .fold(0.0, f64::max)
        .max(
            (0..3)
                .map(|a| spacing[a] * control_box.upper[a] as f64)
                .fold(0.0, f64::max),
        );
    // Conservatively refuse roundoff-near surface points. This is a
    // numerical domain guard, not added physical surface thickness.
    let guard = Dd::from_f64(64.0 * f64::EPSILON * scale);
    let dot = |a: [Dd; 3], b: [Dd; 3]| (0..3).fold(Dd::ZERO, |s, i| s + a[i] * b[i]);
    for axis in 0..3 {
        let a = columns[(axis + 1) % 3];
        let b = columns[(axis + 2) % 3];
        let normal = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        let determinant = dot(normal, columns[axis]);
        let coordinate = dot(normal, delta) / determinant;
        let low = Dd::from_f64(spacing[axis]) * Dd::from_f64(control_box.lower[axis] as f64);
        let high = Dd::from_f64(spacing[axis]) * Dd::from_f64(control_box.upper[axis] as f64);
        if !coordinate.hi.is_finite() || !coordinate.lo.is_finite() {
            return Err(invalid());
        }
        if coordinate.lt(low - guard) || (high + guard).lt(coordinate) {
            return Ok(());
        }
    }
    Err(invalid())
}

#[cfg(test)]
mod tests;
