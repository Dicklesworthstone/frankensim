//! Original joint volume/force/stress histories, with one physical clock.
//!
//! Recording retains the accepted real samples, not Welch factors. Preparation
//! explicitly selects the existing clamped C2 reconstruction, with authored
//! endpoint rates, initial force/stress integrals and a declared source band.
//! The declaration must fit the original Nyquist interval; it does not certify
//! aliasing, derivative accuracy or a noisy CFD source. Those require source
//! and timestep refinement. No resampling or extrapolation is performed.
//!
//! Molecular observation reuses the complete original Q/F/T propagator. It is
//! explicitly initially resting: zero past and unexcited propagation material
//! are declared, and incompatible initial source data refuse. The retained
//! finite prefix does not establish the later acoustic or material tail.

use super::{
    MolecularSampleClock, MolecularSourceRecord, MolecularTransientBudget,
    MolecularTransientConfig, PreparedMolecularSourceObserver, TransientAcousticSource,
};
use crate::AeroacError;
use crate::curle3d::quadrupole::transient::history::{
    CubicHistory, CubicHistoryEnds, CubicSourceSample, QuadrupoleHistoryBudget,
    QuadrupoleHistoryEnds, stress_sample,
};
use crate::curle3d::quadrupole::transient::molecular::{invalid, nonfinite, poll};
use crate::material::MaterialAeroacError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::relaxation::RelaxationState;
use fs_math::dd::Dd;

/// Authored physical conditions for one finite volume/force record.
/// No missing initial impulse or endpoint derivative is inferred from samples.
#[derive(Clone, Copy, Debug)]
pub struct CompactHistoryEnds {
    pub initial_volume_flow_rate_m3_s2: f64,
    pub final_volume_flow_rate_m3_s2: f64,
    pub initial_force_rate_n_s: [f64; 3],
    pub final_force_rate_n_s: [f64; 3],
    pub initial_force_impulse_n_s: [f64; 3],
}

impl CompactHistoryEnds {
    fn channels(self) -> CubicHistoryEnds<4> {
        CubicHistoryEnds {
            initial_rate: [
                self.initial_volume_flow_rate_m3_s2,
                self.initial_force_rate_n_s[0],
                self.initial_force_rate_n_s[1],
                self.initial_force_rate_n_s[2],
            ],
            final_rate: [
                self.final_volume_flow_rate_m3_s2,
                self.final_force_rate_n_s[0],
                self.final_force_rate_n_s[1],
                self.final_force_rate_n_s[2],
            ],
            // The Q integral is unused by the original radiation operator.
            initial_integral: [
                0.0,
                self.initial_force_impulse_n_s[0],
                self.initial_force_impulse_n_s[1],
                self.initial_force_impulse_n_s[2],
            ],
        }
    }
}

/// Bounds original samples and the SAME reconstruction used for observation.
#[derive(Clone, Copy, Debug)]
pub struct AcousticSourceHistoryBudget {
    pub maximum_compact_sources: usize,
    pub maximum_quadrupole_sources: usize,
    pub maximum_samples_per_source: usize,
    /// Bounds (4*S + 6*V)*N. Capture and C2 preparation use O(this) work;
    /// this is a channel/sample extent, not a measured FLOP count.
    pub maximum_component_samples: usize,
    /// Peak backing arrays for this recorder PLUS one prepared history,
    /// including original timestamps, geometry, C2 nodes and preparation work.
    /// Native CFD samplers, propagation/outputs, caller endpoint data, object
    /// metadata and allocator overhead have separate owners/budgets.
    pub maximum_workspace_bytes: usize,
}

/// Bounded sample-major real Q/F and T records with immutable source identity.
///
/// Both groups are nonempty and remain disjoint physical mechanisms. Manual
/// callers own the common-state/source-validity obligation; the native FW-H
/// adapter enforces one accepted flow borrow and exterior-only volume cells.
/// Pushes allocate nothing; a refused/cancelled sample leaves the accepted
/// prefix and clock unchanged. Only the unpublished next row can be dirty.
pub struct AcousticSourceHistoryRecorder {
    acoustic: FluidAcousticState,
    compact_positions_m: Vec<[f64; 3]>,
    stress_positions_m: Vec<[f64; 3]>,
    clock: MolecularSampleClock,
    times_s: Vec<f64>,
    compact: Vec<[f64; 4]>,
    stress: Vec<[f64; 6]>,
    accepted: usize,
    workspace_bytes: usize,
}

impl AcousticSourceHistoryRecorder {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        acoustic: FluidAcousticState,
        compact_positions_m: &[[f64; 3]],
        stress_positions_m: &[[f64; 3]],
        clock: MolecularSampleClock,
        samples: usize,
        budget: AcousticSourceHistoryBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        let (s, v) = (compact_positions_m.len(), stress_positions_m.len());
        if s == 0
            || v == 0
            || s > budget.maximum_compact_sources
            || v > budget.maximum_quadrupole_sources
            || samples < 2
            || samples > budget.maximum_samples_per_source
        {
            return Err(invalid("joint source history count or sample budget"));
        }
        clock.validate_last_index(samples - 1)?;
        let channels = s
            .checked_mul(4)
            .and_then(|n| n.checked_add(v.checked_mul(6)?))
            .ok_or_else(|| invalid("joint source history channel extent overflow"))?;
        let component_samples = channels
            .checked_mul(samples)
            .filter(|&n| n <= budget.maximum_component_samples)
            .ok_or_else(|| invalid("joint source history work extent budget"))?;
        let geometry = s
            .checked_add(v)
            .and_then(|n| n.checked_mul(size_of::<[f64; 3]>()))
            .ok_or_else(|| invalid("joint source history geometry extent overflow"))?;
        let compact_workspace = CubicHistory::<4>::required_workspace_bytes(s, samples)?;
        let stress_workspace = CubicHistory::<6>::required_workspace_bytes(v, samples)?;
        let workspace_bytes = component_samples
            .checked_add(samples)
            .and_then(|n| n.checked_mul(size_of::<f64>()))
            // Recorder positions and prepared combined source positions. Each
            // original C2 owner also accounts for its retained position slice.
            .and_then(|n| n.checked_add(geometry.checked_mul(2)?))
            .and_then(|n| n.checked_add(compact_workspace))
            .and_then(|n| n.checked_add(stress_workspace))
            .filter(|&n| n <= budget.maximum_workspace_bytes)
            .ok_or_else(|| invalid("joint source history storage budget"))?;
        let compact_positions_m = copy_positions(compact_positions_m, &mut cancelled)?;
        let stress_positions_m = copy_positions(stress_positions_m, &mut cancelled)?;
        // Both products are bounded by the admitted channel/sample extent.
        let times_s = allocate(samples, 0.0, &mut cancelled)?;
        let compact = allocate(samples * s, [0.0; 4], &mut cancelled)?;
        let stress = allocate(samples * v, [0.0; 6], &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(Self {
            acoustic,
            compact_positions_m,
            stress_positions_m,
            clock,
            times_s,
            compact,
            stress,
            accepted: 0,
            workspace_bytes,
        })
    }

    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }
    pub const fn clock(&self) -> MolecularSampleClock {
        self.clock
    }
    pub fn compact_source_positions_m(&self) -> &[[f64; 3]] {
        &self.compact_positions_m
    }
    pub fn quadrupole_source_positions_m(&self) -> &[[f64; 3]] {
        &self.stress_positions_m
    }
    pub const fn accepted_sample_count(&self) -> usize {
        self.accepted
    }
    pub fn remaining_sample_count(&self) -> usize {
        self.times_s.len() - self.accepted
    }
    pub fn sample_times_s(&self) -> &[f64] {
        &self.times_s[..self.accepted]
    }
    /// Accepted sample-major `[Q,Fx,Fy,Fz]`, with original SI source scaling.
    pub fn volume_force_samples(&self) -> &[[f64; 4]] {
        &self.compact[..self.accepted * self.compact_positions_m.len()]
    }
    /// Accepted sample-major `[xx,yy,zz,xy,xz,yz]`; shears are not doubled.
    pub fn stress_samples_n_m(&self) -> &[[f64; 6]] {
        &self.stress[..self.accepted * self.stress_positions_m.len()]
    }
    /// Admitted peak including one simultaneous immutable reconstruction.
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    #[allow(clippy::too_many_arguments)]
    pub fn push_sample(
        &mut self,
        compact_positions_m: &[[f64; 3]],
        stress_positions_m: &[[f64; 3]],
        time_s: f64,
        volume_force: &[[f64; 4]],
        stress: &[[f64; 6]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        poll(&mut cancelled)?;
        if self.remaining_sample_count() == 0
            || volume_force.len() != self.compact_positions_m.len()
            || stress.len() != self.stress_positions_m.len()
            || compact_positions_m.len() != self.compact_positions_m.len()
            || stress_positions_m.len() != self.stress_positions_m.len()
        {
            return Err(invalid("joint source history sample shape or full record"));
        }
        self.clock.admit(self.accepted, time_s)?;
        for (i, (given, original)) in compact_positions_m
            .iter()
            .chain(stress_positions_m)
            .zip(
                self.compact_positions_m
                    .iter()
                    .chain(&self.stress_positions_m),
            )
            .enumerate()
        {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            if given != original {
                return Err(invalid("joint source history geometry/order changed"));
            }
        }
        let start = self.accepted * volume_force.len();
        for (i, row) in volume_force.iter().enumerate() {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            if row.iter().any(|x| !x.is_finite()) {
                return Err(nonfinite("joint source history volume/force sample"));
            }
            self.compact[start + i] = *row;
        }
        let start = self.accepted * stress.len();
        for (i, row) in stress.iter().enumerate() {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            if row.iter().any(|x| !x.is_finite()) {
                return Err(nonfinite("joint source history stress sample"));
            }
            self.stress[start + i] = *row;
        }
        poll(&mut cancelled)?;
        self.times_s[self.accepted] = time_s;
        self.accepted += 1;
        Ok(())
    }

    /// Select one offline clamped C2 model for ALL original physical channels.
    ///
    /// `maximum_source_frequency_hz` is an explicit caller declaration, finite
    /// and strictly below the ORIGINAL Nyquist frequency. It is neither a
    /// filter nor evidence that the recorded source satisfies that bandwidth.
    /// Endpoint conditions and cumulative integrals remain physical inputs.
    /// This borrows the completed record: refusal/cancellation is retryable.
    pub fn prepare_clamped_cubic(
        &self,
        compact_ends: &[CompactHistoryEnds],
        stress_ends: &[QuadrupoleHistoryEnds],
        maximum_source_frequency_hz: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAcousticSourceHistory, AeroacError> {
        poll(&mut cancelled)?;
        if self.remaining_sample_count() != 0
            || compact_ends.len() != self.compact_positions_m.len()
            || stress_ends.len() != self.stress_positions_m.len()
        {
            return Err(invalid(
                "joint source history incomplete record or endpoint shape",
            ));
        }
        let normalized_band = maximum_source_frequency_hz * self.clock.timestep_s();
        if !maximum_source_frequency_hz.is_finite()
            || maximum_source_frequency_hz <= 0.0
            || !normalized_band.is_finite()
            || normalized_band <= 0.0
            || normalized_band >= 0.5
        {
            return Err(invalid(
                "declared source bandwidth must be below original Nyquist",
            ));
        }
        let s = self.compact_positions_m.len();
        let v = self.stress_positions_m.len();
        let make_budget = |sources, bytes| QuadrupoleHistoryBudget {
            maximum_sources: sources,
            maximum_samples_per_source: self.accepted,
            maximum_workspace_bytes: bytes,
        };
        let compact = CubicHistory::new(
            &self.compact_positions_m,
            self.clock,
            &self.times_s,
            |sample, source| self.compact[sample * s + source],
            |source| compact_ends[source].channels(),
            make_budget(
                s,
                CubicHistory::<4>::required_workspace_bytes(s, self.accepted)?,
            ),
            &mut cancelled,
        )?;
        let stress = CubicHistory::new(
            &self.stress_positions_m,
            self.clock,
            &self.times_s,
            |sample, source| self.stress[sample * v + source],
            |source| CubicHistoryEnds {
                initial_rate: stress_ends[source].initial_rate_n_m_s,
                final_rate: stress_ends[source].final_rate_n_m_s,
                initial_integral: stress_ends[source].initial_integral_n_m_s,
            },
            make_budget(
                v,
                CubicHistory::<6>::required_workspace_bytes(v, self.accepted)?,
            ),
            &mut cancelled,
        )?;
        let mut positions = reserve(s + v)?;
        for (i, &x) in self
            .compact_positions_m
            .iter()
            .chain(&self.stress_positions_m)
            .enumerate()
        {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            positions.push(x);
        }
        let workspace_bytes = compact.workspace_bytes()
            + stress.workspace_bytes()
            + positions.len() * size_of::<[f64; 3]>();
        poll(&mut cancelled)?;
        Ok(PreparedAcousticSourceHistory {
            acoustic: self.acoustic,
            positions,
            compact,
            stress,
            maximum_source_frequency_hz,
            workspace_bytes,
        })
    }
}

/// Immutable original-layout history with consistent derivatives and integrals.
/// Source indices are compact positions first, then exterior stress positions;
/// absent mechanisms are exactly zero, with no merging of coincident sources.
pub struct PreparedAcousticSourceHistory {
    acoustic: FluidAcousticState,
    positions: Vec<[f64; 3]>,
    compact: CubicHistory<4>,
    stress: CubicHistory<6>,
    maximum_source_frequency_hz: f64,
    workspace_bytes: usize,
}

impl PreparedAcousticSourceHistory {
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }
    pub const fn clock(&self) -> MolecularSampleClock {
        self.compact.clock()
    }
    pub fn source_positions_m(&self) -> &[[f64; 3]] {
        &self.positions
    }
    pub fn compact_source_count(&self) -> usize {
        self.compact.source_count()
    }
    pub fn quadrupole_source_count(&self) -> usize {
        self.stress.source_count()
    }
    pub const fn sample_count(&self) -> usize {
        self.compact.sample_count()
    }
    pub fn source_time_bounds_s(&self) -> [Dd; 2] {
        self.compact.source_time_bounds_s()
    }
    pub const fn declared_maximum_source_frequency_hz(&self) -> f64 {
        self.maximum_source_frequency_hz
    }
    /// Requested reconstruction payload including temporary preparation work.
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    pub fn initial_source(&self, source: usize) -> Result<TransientAcousticSource, AeroacError> {
        if source < self.compact.source_count() {
            self.compact.initial_source(source).map(compact_sample)
        } else {
            self.stress
                .initial_source(source - self.compact.source_count())
                .map(|sample| TransientAcousticSource {
                    stress: stress_sample(sample),
                    ..TransientAcousticSource::default()
                })
        }
    }

    /// Query the original finite history, keeping retarded times in DD seconds.
    pub fn source_at_clock_dd(
        &self,
        source: usize,
        time_s: Dd,
    ) -> Result<TransientAcousticSource, AeroacError> {
        if source < self.compact.source_count() {
            self.compact
                .source_at_clock_dd(source, time_s)
                .map(compact_sample)
        } else {
            self.stress
                .source_at_clock_dd(source - self.compact.source_count(), time_s)
                .map(|sample| TransientAcousticSource {
                    stress: stress_sample(sample),
                    ..TransientAcousticSource::default()
                })
        }
    }

    /// Bind the original equilibrium reference to the SAME material's causal
    /// source propagator. Its frozen-speed front stays inside the native owner.
    pub fn prepare_molecular_observer(
        &self,
        molecular: RelaxationState,
        receiver_m: [f64; 3],
        config: MolecularTransientConfig,
        budget: MolecularTransientBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularSourceObserver, MaterialAeroacError> {
        poll(&mut cancelled)?;
        if molecular
            .gas()
            .acoustic_state()
            .map_err(MaterialAeroacError::Material)?
            != self.acoustic
            || config.timestep_s != self.clock().timestep_s()
            || config.samples > self.sample_count()
        {
            return Err(
                invalid("molecular source history reference, clock or horizon mismatch").into(),
            );
        }
        PreparedMolecularSourceObserver::new(
            molecular,
            receiver_m,
            &self.positions,
            config,
            budget,
            cancelled,
        )
    }

    /// Declare zero past and an initially unexcited exterior/material, then
    /// propagate the complete coherent field with the original CQ owner.
    /// Original values, force impulse and stress value/rate/integral must have
    /// the compatible zero onset required by that owner. No demeaning, onset
    /// repair, missing impulse, tail extension or second propagation delay.
    /// Exact matching ideal-image/wedge paths may be supplied by `from_paths`.
    pub fn observe_initially_resting(
        &self,
        observer: &PreparedMolecularSourceObserver,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularSourceRecord, AeroacError> {
        poll(&mut cancelled)?;
        if observer.paths().source_count() != self.positions.len()
            || observer.config().samples > self.sample_count()
            || observer.config().timestep_s != self.clock().timestep_s()
            || !observer
                .molecular_state()
                .gas()
                .acoustic_state()
                .is_ok_and(|state| state == self.acoustic)
        {
            return Err(invalid(
                "molecular observer changed original history reference or extent",
            ));
        }
        for (i, &x) in self.positions.iter().enumerate() {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            if observer.paths().source_position_m(i) != Some(x) {
                return Err(invalid(
                    "molecular observer changed original history geometry/order",
                ));
            }
        }
        let epoch = Dd::from_f64(self.clock().epoch_s());
        observer.observe_initially_resting(
            self.clock(),
            |source, time| {
                if time == epoch {
                    self.initial_source(source)
                } else {
                    self.source_at_clock_dd(source, time)
                }
            },
            cancelled,
        )
    }
}

fn compact_sample(sample: CubicSourceSample<4>) -> TransientAcousticSource {
    let mut source = TransientAcousticSource::default();
    source.compact.volume_flow_m3_s = sample.value[0];
    source.compact.volume_flow_rate_m3_s2 = sample.rate[0];
    for i in 0..3 {
        source.compact.force_on_fluid_n[i] = sample.value[i + 1];
        source.compact.force_rate_on_fluid_n_s[i] = sample.rate[i + 1];
        source.force_impulse_n_s[i] = sample.integral[i + 1];
    }
    source
}

fn reserve<T>(count: usize) -> Result<Vec<T>, AeroacError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| invalid("joint source history allocation exceeds available memory"))?;
    Ok(values)
}

fn allocate<T: Copy>(
    count: usize,
    value: T,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<T>, AeroacError> {
    let mut values = reserve(count)?;
    for i in 0..count {
        if i.is_multiple_of(256) {
            poll(cancelled)?;
        }
        values.push(value);
    }
    Ok(values)
}

fn copy_positions(
    positions: &[[f64; 3]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<[f64; 3]>, AeroacError> {
    let mut copy = reserve(positions.len())?;
    for (i, &position) in positions.iter().enumerate() {
        if i.is_multiple_of(256) {
            poll(cancelled)?;
        }
        if position.iter().any(|x| !x.is_finite()) {
            return Err(nonfinite("joint source history geometry"));
        }
        copy.push(position);
    }
    Ok(copy)
}

#[cfg(test)]
mod tests;
