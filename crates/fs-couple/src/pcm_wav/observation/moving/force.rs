//! Causal full acoustic fields of accepted moving force and volume sources.
//!
//! Every path integrates the original force at the original moving geometry:
//! J = integral (3 n (n.F) - F)/R^3 ds. Immutable prefixes split at BOTH accepted
//! force knots and accepted trajectory knots. Fixed microphones cache the
//! complete cells and evaluate only the remaining partial cell. Moving
//! microphones integrate the SAME partition at their current spatial point.
//! Local retarded terms come from fs-aeroac's moving-force field owner.
//!
//! This owner declares zero force before the first strength knot and requires
//! zero force AND force rate there. It adds no unknown initial acoustic field.
//! Volume-source prehistory retains the existing local-potential contract.
//! Source bandwidth, compactness, interpolation convergence and linearity remain
//! producer obligations. Received bandwidth must include geometric modulation.
//! The quadrature target is an adaptive error estimate, not an interval proof.
//! A continuous-speed geometric refinement guard additionally prevents a
//! narrow close pass from falling unnoticed between quadrature samples.

pub(in crate::pcm_wav::observation) mod integral;
mod motion;
mod receiver;

pub use integral::ForceHistoryQuadrature;
pub use motion::{FixedForceMotion, ForceHistoryPath, ForceMotion};

use super::field::CompactFieldSample;
use super::*;
use fs_bem::near_field::FirstOrder;
use fs_material::fluid_acoustics::FluidAcousticState;
use integral::Work;

#[derive(Clone, Copy)]
struct HistoryPrefix {
    source_time_s: f64,
    strength_interval: usize,
    integral: [Dd; 3],
}

struct ForcePath {
    source_index: usize,
    path_index: usize,
    geometry: ForceHistoryPath,
    prefix_start: usize,
    prefix_count: usize,
}

/// Full moving-dipole/monopole fields at fixed or moving microphones.
///
/// Uniform and accepted curved sources share the original pressure histories,
/// retarded clocks and Doppler admission. Ideal-wall paths each retain their
/// own independently retarded image path. Pressure and transformed vector velocity
/// are summed coherently before intensity or directional microphone pickup.
///
/// Construction charges prefixes, path records and field scratch to the
/// original payload budget. Quadrature work is bounded across preparation AND
/// callback attempts. Callbacks allocate nothing. A refused or cancelled block
/// publishes neither samples nor clock, but actual quadrature work is not
/// refunded. Successful output is independent of block size.
///
/// Native moving source/receiver pairs reuse the original motion and history
/// owners, including tangential uniform flow and complete ideal image groups.
/// Their spatial history kernel changes with each microphone event, so the
/// entire past-force integral is recomputed under the same lifetime work cap.
/// Continuous source-speed bounds enforce past-worldline separation during
/// integration; a newly encountered singular past refuses the entire block.
/// Homogeneous initial acoustic fields compose separately.
pub struct MovingForceRenderer<'a, M: ForceMotion = UniformSourceMotion> {
    original: MovingCompactRenderer<'a, M>,
    paths: Vec<ForcePath>,
    prefixes: Vec<HistoryPrefix>,
    scratch: Vec<CompactFieldSample>,
    microphone: FirstOrder,
    medium: FluidAcousticState,
    quadrature: ForceHistoryQuadrature,
    work: Work,
}

impl<'a, M: ForceMotion> MovingForceRenderer<'a, M> {
    /// Admit original past source coverage and explicit zero force prehistory.
    ///
    /// A nonzero initial force would imply an unrepresented onset wavefront.
    /// Requiring both initial force and derivative zero also preserves C1 onset.
    /// Fixed receivers admit complete past geometry now; moving receivers check
    /// past-worldline clearance at each current image point before publishing
    /// an output block.
    pub fn new_with_zero_force_prehistory(
        observer: M::Observer,
        sources: &[CompactSourceHistory<'a, M>],
        microphone: FirstOrder,
        config: MovingPressureConfig,
        mut budget: MovingPressureBudget,
        quadrature: ForceHistoryQuadrature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        let mut work = Work::new(quadrature)?;
        let field_bytes = config
            .max_block
            .checked_mul(size_of::<CompactFieldSample>())
            .ok_or_else(|| sizing("moving force field scratch payload overflows"))?;
        budget.maximum_payload_bytes = budget
            .maximum_payload_bytes
            .checked_sub(field_bytes)
            .ok_or_else(|| sizing("moving force fields exceed the retained payload budget"))?;
        let original = MovingCompactRenderer::new(observer, sources, config, budget, cancelled)?;
        for track in &original.tracks {
            let first = track.source.samples[0].strength;
            if first.force_on_fluid_n != [0.0; 3] || first.force_rate_on_fluid_n_s != [0.0; 3] {
                return Err(sizing(
                    "zero force prehistory requires zero initial force and force rate",
                ));
            }
        }
        let physical_observer = original.observer_at(config.first_receiver_time_s)?;
        let medium = physical_observer.medium();
        let original_bytes = original
            .tracks
            .len()
            .checked_mul(size_of::<PreparedTrack<'a, M>>())
            .and_then(|n| {
                original
                    .intervals
                    .len()
                    .checked_mul(size_of::<CubicInterval>())
                    .and_then(|m| n.checked_add(m))
            })
            .and_then(|n| {
                original
                    .scratch
                    .len()
                    .checked_mul(size_of::<f64>())
                    .and_then(|m| n.checked_add(m))
            })
            .ok_or_else(|| sizing("moving force base payload overflows"))?;
        let path_bytes = original
            .source_path_count()
            .checked_mul(size_of::<ForcePath>())
            .ok_or_else(|| sizing("moving force path payload overflows"))?;
        let prefix_bytes = budget
            .maximum_payload_bytes
            .checked_sub(original_bytes)
            .and_then(|n| n.checked_sub(path_bytes))
            .ok_or_else(|| sizing("moving force paths exceed the retained payload budget"))?;
        let maximum_prefixes = prefix_bytes / size_of::<HistoryPrefix>();
        let mut paths = Vec::new();
        paths
            .try_reserve_exact(original.source_path_count())
            .map_err(|_| sizing("cannot reserve moving force paths"))?;
        let last_time = original.receiver_time_at(config.samples - 1)?;
        let mut prefix_count = 0_usize;
        for (source_index, track) in original.tracks.iter().enumerate() {
            check_cancel(cancelled)?;
            for path_index in 0..track.source.motion.path_count() {
                let last = prepare_path_history(
                    &original.observer,
                    &track.source,
                    path_index,
                    last_time,
                    config.maximum_time_residual_s,
                    cancelled,
                )?;
                let first_time = track.source.samples[0].source_time_s;
                let last_strength =
                    track.source.samples[track.source.samples.len() - 1].source_time_s;
                // Arrival slope >= 1-Mmax. Both the final and a queried root
                // may use the residual allowance. Round every bound outward.
                let source_observer = track
                    .source
                    .motion
                    .observer_at(&original.observer, last_time)?;
                let denominator = (1.0 - source_observer.limits().maximum_mach).next_down();
                let padding =
                    ((2.0 * config.maximum_time_residual_s).next_up() / denominator).next_up();
                let covered_end = (last.source_time_s() + padding)
                    .next_up()
                    .min(last_strength);
                let geometry = track.source.motion.force_history_path_at(
                    &original.observer,
                    path_index,
                    config.first_receiver_time_s,
                    [first_time, covered_end],
                    cancelled,
                )?;
                let duration = geometry.source_end_s - first_time;
                if !duration.is_finite() || duration < 0.0 {
                    return Err(sizing("complete force history duration is unrepresentable"));
                }
                let prefix_start = prefix_count;
                prefix_count = charge_prefix(prefix_count, maximum_prefixes)?;
                let mut time = first_time;
                let mut strength_interval = 0;
                while let Some((next, next_interval)) = next_breakpoint(
                    &track.source,
                    time,
                    strength_interval,
                    geometry.source_end_s,
                )? {
                    check_cancel(cancelled)?;
                    prefix_count = charge_prefix(prefix_count, maximum_prefixes)?;
                    time = next;
                    strength_interval = next_interval;
                }
                paths.push(ForcePath {
                    source_index,
                    path_index,
                    geometry,
                    prefix_start,
                    prefix_count: prefix_count - prefix_start,
                });
            }
        }
        let mut prefixes = Vec::new();
        prefixes
            .try_reserve_exact(prefix_count)
            .map_err(|_| sizing("cannot reserve moving force history prefixes"))?;
        for path in &paths {
            check_cancel(cancelled)?;
            let track = &original.tracks[path.source_index];
            let first = track.source.samples[0].source_time_s;
            let duration = path.geometry.source_end_s - first;
            let mut current = HistoryPrefix {
                source_time_s: first,
                strength_interval: 0,
                integral: [Dd::ZERO; 3],
            };
            prefixes.push(current);
            while let Some((next, next_interval)) = next_breakpoint(
                &track.source,
                current.source_time_s,
                current.strength_interval,
                path.geometry.source_end_s,
            )? {
                let increment = if M::MOVING_HISTORY_KERNEL {
                    // Store only immutable knot metadata. A prefix evaluated
                    // at this microphone event is invalid at the next event.
                    [Dd::ZERO; 3]
                } else {
                    integrate_cell(
                        track,
                        &original.intervals,
                        path.geometry,
                        current.strength_interval,
                        [current.source_time_s, next],
                        duration,
                        quadrature,
                        &mut work,
                        cancelled,
                    )?
                };
                current = HistoryPrefix {
                    source_time_s: next,
                    strength_interval: next_interval,
                    integral: add_integrals(current.integral, increment)?,
                };
                prefixes.push(current);
            }
        }
        let mut scratch = Vec::new();
        scratch
            .try_reserve_exact(config.max_block)
            .map_err(|_| sizing("cannot reserve moving force field scratch"))?;
        scratch.resize(config.max_block, CompactFieldSample::default());
        check_cancel(cancelled)?;
        Ok(Self {
            original,
            paths,
            prefixes,
            scratch,
            microphone,
            medium,
            quadrature,
            work,
        })
    }

    /// Original physical receiver and boundary identity.
    pub const fn observer(&self) -> &M::Observer {
        self.original.observer()
    }
    /// Original observation clock, bandwidth and finite horizon.
    pub const fn config(&self) -> MovingPressureConfig {
        self.original.config()
    }
    /// Fixed physical first-order microphone pickup.
    pub const fn microphone(&self) -> FirstOrder {
        self.microphone
    }

    /// Background fluid velocity in the original laboratory frame [m/s].
    pub fn mean_flow_velocity_m_s(&self) -> [f64; 3] {
        self.original.mean_flow_velocity_m_s()
    }

    /// Energy and laboratory transport after coherent source/image summation.
    /// The original fluid/frame supplies rho, c and U; microphone pickup and
    /// microphone translation are not acoustic fields or fluid advection.
    pub fn energy_transport(
        &self,
        field: &CompactFieldSample,
    ) -> Result<fs_aeroac::curle3d::energy::AcousticEnergyTransport, RenderError> {
        fs_aeroac::curle3d::energy::instantaneous(
            self.medium,
            self.mean_flow_velocity_m_s(),
            field.pressure_pa,
            field.particle_velocity_m_s,
        )
        .map_err(acoustic)
    }
    /// Original physical source count.
    pub fn source_count(&self) -> usize {
        self.original.source_count()
    }
    /// All original direct and image paths.
    pub const fn source_path_count(&self) -> usize {
        self.original.source_path_count()
    }
    /// Complete samples published on the shared receiver clock.
    pub const fn samples_rendered(&self) -> u64 {
        self.original.samples_rendered()
    }
    /// Remaining original receiver samples.
    pub fn remaining_samples(&self) -> u64 {
        self.original.remaining_samples()
    }
    /// Largest callback fitting retained scratch.
    pub const fn max_block_len(&self) -> usize {
        self.original.max_block_len()
    }
    /// Actual integrand evaluations, including failed/cancelled attempts.
    pub const fn quadrature_evaluations(&self) -> u64 {
        self.work.used
    }
    /// Original integer-index timestamp, independent of block partition.
    pub fn receiver_time_at(&self, sample: u64) -> Result<f64, RenderError> {
        self.original.receiver_time_at(sample)
    }

    /// Publish coherent pressure, velocity, intensity and pickup atomically.
    pub fn field_block(&mut self, output: &mut [CompactFieldSample]) -> Result<(), RenderError> {
        self.field_block_controlled(output, &mut || false)
    }
    /// Full-field callback with bounded cancellation polling.
    pub fn field_block_controlled(
        &mut self,
        output: &mut [CompactFieldSample],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        self.prepare_block(output.len(), cancelled)?;
        output.copy_from_slice(&self.scratch[..output.len()]);
        self.original.completed += output.len() as u64;
        Ok(())
    }
    /// Emit Pa-equivalent directional microphone output.
    pub fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.block_controlled(output, &mut || false)
    }
    /// Atomic audio callback on the same sample cursor as full-field output.
    pub fn block_controlled(
        &mut self,
        output: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        self.prepare_block(output.len(), cancelled)?;
        for (out, field) in output.iter_mut().zip(&self.scratch) {
            *out = field.microphone_pressure_pa;
        }
        self.original.completed += output.len() as u64;
        Ok(())
    }

    fn prepare_block(
        &mut self,
        count: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        if count == 0 {
            return Err(RenderError::EmptyBlock);
        }
        if count > self.max_block_len() {
            return Err(sizing("moving force field block exceeds its capacity"));
        }
        self.validate_sample_count(
            u64::try_from(count).map_err(|_| sizing("moving force field count exceeds u64"))?,
        )?;
        for i in 0..count {
            check_cancel(cancelled)?;
            let time = self.receiver_time_at(self.original.completed + i as u64)?;
            let mut pressure = Dd::ZERO;
            let mut velocity = [Dd::ZERO; 3];
            for path in &self.paths {
                check_cancel(cancelled)?;
                let track = &self.original.tracks[path.source_index];
                let prepared = prepare_path_history(
                    &self.original.observer,
                    &track.source,
                    path.path_index,
                    time,
                    self.original.config.maximum_time_residual_s,
                    cancelled,
                )?;
                let source_time = prepared.source_time_s();
                let bank = &self.prefixes[path.prefix_start..path.prefix_start + path.prefix_count];
                let upper = bank.partition_point(|prefix| prefix.source_time_s <= source_time);
                if upper == 0 || source_time > path.geometry.source_end_s {
                    return Err(sizing(
                        "retarded force query exceeds its admitted past geometry",
                    ));
                }
                let first_time = track.source.samples[0].source_time_s;
                let geometry = if M::MOVING_HISTORY_KERNEL {
                    track.source.motion.force_history_path_at(
                        &self.original.observer,
                        path.path_index,
                        time,
                        [first_time, path.geometry.source_end_s],
                        cancelled,
                    )?
                } else {
                    path.geometry
                };
                let duration = geometry.source_end_s - first_time;
                let prefix = bank[upper - 1];
                let mut integral = prefix.integral;
                if M::MOVING_HISTORY_KERNEL {
                    // Every past cell sees the SAME current virtual receiver.
                    // Original source/force knots are shared across output
                    // events; no interpolation or trajectory is reconstructed.
                    for cell in bank[..upper].windows(2) {
                        let increment = integrate_cell(
                            track,
                            &self.original.intervals,
                            geometry,
                            cell[0].strength_interval,
                            [cell[0].source_time_s, cell[1].source_time_s],
                            duration,
                            self.quadrature,
                            &mut self.work,
                            cancelled,
                        )?;
                        integral = add_integrals(integral, increment)?;
                    }
                }
                let increment = integrate_cell(
                    track,
                    &self.original.intervals,
                    geometry,
                    prefix.strength_interval,
                    [prefix.source_time_s, source_time],
                    duration,
                    self.quadrature,
                    &mut self.work,
                    cancelled,
                )?;
                let integral = add_integrals(integral, increment)?;
                let strength = interpolate_interval(
                    track,
                    &self.original.intervals,
                    prefix.strength_interval,
                    source_time,
                )?;
                let kinematics = track
                    .source
                    .motion
                    .kinematics_at(source_time)
                    .map_err(acoustic)?;
                let (_, path_pressure, path_velocity) = geometry
                    .observer
                    .prepare_force_emission(source_time, kinematics)
                    .and_then(|local| local.evaluate_dd(strength, integral))
                    .map_err(acoustic)?;
                pressure = pressure + Dd::from_f64(geometry.pressure_sign) * path_pressure;
                for axis in 0..3 {
                    let component = match geometry.velocity_transform {
                        None => path_velocity[axis],
                        Some(matrix) => (0..3).fold(Dd::ZERO, |sum, j| {
                            sum + Dd::from_f64(matrix[axis][j]) * path_velocity[j]
                        }),
                    };
                    velocity[axis] = velocity[axis] + component;
                }
            }
            self.scratch[i] = CompactFieldSample::observe(
                time,
                pressure,
                velocity,
                self.microphone,
                self.medium,
            )?;
        }
        check_cancel(cancelled)
    }
}

impl<M: ForceMotion> PressureRenderer for MovingForceRenderer<'_, M> {
    fn samples_rendered(&self) -> u64 {
        Self::samples_rendered(self)
    }
    fn max_block_len(&self) -> usize {
        Self::max_block_len(self)
    }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.original.validate_sample_rate(rate)
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        self.original.validate_sample_count(samples)
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        Self::block(self, output)
    }
}

fn charge_prefix(count: usize, maximum: usize) -> Result<usize, RenderError> {
    count
        .checked_add(1)
        .filter(|next| *next <= maximum)
        .ok_or_else(|| sizing("moving force prefixes exceed the retained payload budget"))
}

// There is deliberately no prefix at an artificial padded horizon: only real
// strength/geometry knots split the interpolants. Runtime integrates the last
// partial cell, avoiding tiny artificial cells and preserving one-cell work.
fn next_breakpoint<M: ForceMotion>(
    source: &CompactSourceHistory<'_, M>,
    time: f64,
    strength_interval: usize,
    end: f64,
) -> Result<Option<(f64, usize)>, RenderError> {
    if time >= end {
        return Ok(None);
    }
    let next_strength = source
        .samples
        .get(strength_interval + 1)
        .map(|sample| sample.source_time_s)
        .unwrap_or(f64::INFINITY);
    let next_geometry = source
        .motion
        .next_force_history_knot_after(time)?
        .unwrap_or(f64::INFINITY);
    let next = next_strength.min(next_geometry);
    if next > end {
        return Ok(None);
    }
    if !next.is_finite() || next <= time {
        return Err(sizing(
            "force prefix partition lost its original source time order",
        ));
    }
    Ok(Some((
        next,
        strength_interval + usize::from(next == next_strength),
    )))
}

fn add_integrals(first: [Dd; 3], second: [Dd; 3]) -> Result<[Dd; 3], RenderError> {
    let total: [Dd; 3] = core::array::from_fn(|i| first[i] + second[i]);
    if total.iter().any(|value| !value.to_f64().is_finite()) {
        return Err(sizing("moving force history prefix is unrepresentable"));
    }
    Ok(total)
}

fn integrate_cell<M: ForceMotion>(
    track: &PreparedTrack<'_, M>,
    intervals: &[CubicInterval],
    geometry: ForceHistoryPath,
    strength_interval: usize,
    source_interval: [f64; 2],
    history_duration_s: f64,
    quadrature: ForceHistoryQuadrature,
    work: &mut Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Dd; 3], RenderError> {
    if source_interval[0] == source_interval[1] {
        check_cancel(cancelled)?;
        return Ok([Dd::ZERO; 3]);
    }
    let tolerance = quadrature.absolute_tolerance_n_s_m3
        * ((source_interval[1] - source_interval[0]) / history_duration_s);
    integral::integrate(
        source_interval[0],
        source_interval[1],
        tolerance,
        geometry.maximum_source_speed_m_s,
        geometry.minimum_history_distance_m,
        quadrature,
        work,
        |time| {
            let strength = interpolate_interval(track, intervals, strength_interval, time)?;
            let source = track.source.motion.kinematics_at(time).map_err(acoustic)?;
            let displacement = core::array::from_fn(|axis| {
                geometry.observer.observer_position_m()[axis] - source.position_m[axis]
            });
            let value = geometry
                .observer
                .force_history_integrand(source.position_m, strength.force_on_fluid_n)
                .map_err(acoustic)?;
            Ok((value, norm(displacement)))
        },
        cancelled,
    )
}
