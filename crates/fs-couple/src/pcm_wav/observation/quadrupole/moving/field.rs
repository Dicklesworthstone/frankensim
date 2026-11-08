//! Full moving stress fields from one original finite stress record.
//!
//! Pressure and local velocity use the same retarded event. The additional
//! velocity history is integrated against the ORIGINAL moving source position
//! and current virtual receiver, with the existing bounded vector quadrature.
//! Fixed receivers cache complete source/geometry cells; moving receivers
//! reevaluate those cells at their current point. No plane-wave p/Z surrogate,
//! independent stress reconstruction or inferred prehistory is introduced.

use super::super::{QuadrupoleFieldSample, QuadrupolePressureBudget, check_cancel, sizing};
use super::{MovingQuadrupoleRenderer, MovingQuadrupoleSourceWindow, QuadrupoleMotion};
use crate::pcm_wav::observation::moving::force::{
    ForceHistoryPath, ForceHistoryQuadrature, ForceMotion,
    integral::{self, Work},
};
use crate::pcm_wav::observation::{PressureRenderer, moving::MovingPressureConfig};
use crate::render::RenderError;
use core::mem::size_of;
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::UniformSourceMotion;
use fs_aeroac::curle3d::quadrupole::transient::{
    TransientQuadrupoleSource, history::PreparedQuadrupoleHistory,
};
use fs_bem::near_field::FirstOrder;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{dd::Dd, det};

/// Original pressure capacities plus the retained causal integration payload.
#[derive(Clone, Copy, Debug)]
pub struct MovingQuadrupoleFieldBudget {
    /// The scratch limit covers full field samples. Source histories are borrowed.
    pub pressure: QuadrupolePressureBudget,
    /// Path records and all cached knot prefixes, excluding allocator metadata.
    pub maximum_history_bytes: usize,
    /// Shared three-vector quadrature: the stress kernel also has units N s/m³.
    /// Work includes preparation and failed callback attempts; it is never refunded.
    pub quadrature: ForceHistoryQuadrature,
}

#[derive(Clone, Copy)]
struct Prefix {
    time: f64,
    source_knot: usize,
    integral: [Dd; 3],
}

struct Path {
    source: usize,
    image: usize,
    geometry: ForceHistoryPath,
    start: usize,
    count: usize,
}

/// Coherent moving-quadrupole pressure, fluid velocity, intensity and pickup.
///
/// The original record supplies consistent Q/Qdot/Qddot. Its authored initial
/// Q, Qdot and integral must all vanish: this constructor explicitly declares
/// zero stress before the record and supplies no homogeneous initial field.
/// Nonzero initial Qddot is allowed, with its physical onset wavefront.
/// Every motion starts at its original record position at the record epoch.
///
/// Uniform/curved source motion, moving microphones and ideal images retain
/// the original admission, clocks and geometric limits. The complete field is
/// in the stationary fluid frame. Its pressure-work flux is signed p*u; the
/// fixed-axis first-order pickup is an ideal microphone, not capsule dynamics.
/// Source compactness, bandwidth and spline convergence remain producer duties.
///
/// Callback scratch is retained once. Blocks allocate nothing and publish both
/// output and sample cursor only after success. Actual quadrature work remains
/// charged after refusal or cancellation. The estimate from adaptive quadrature
/// is not an interval error certificate or a source-model validity claim.
pub struct MovingQuadrupoleFieldRenderer<
    'a,
    M: QuadrupoleMotion + ForceMotion = UniformSourceMotion,
> {
    observer: M::Observer,
    sources: &'a [MovingQuadrupoleSourceWindow<M>],
    history: &'a PreparedQuadrupoleHistory,
    config: MovingPressureConfig,
    microphone: FirstOrder,
    medium: FluidAcousticState,
    paths: Vec<Path>,
    prefixes: Vec<Prefix>,
    scratch: Vec<QuadrupoleFieldSample>,
    quadrature: ForceHistoryQuadrature,
    work: Work,
    completed: u64,
    maximum_retarded_iterations: u64,
}

impl<'a, M: QuadrupoleMotion + ForceMotion> MovingQuadrupoleFieldRenderer<'a, M> {
    /// Bind the original finite stress record with explicit zero prehistory.
    /// The supplied windows retain their existing source order and bandwidth;
    /// each is contained in the original record, including all image queries.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_zero_stress_prehistory(
        observer: M::Observer,
        sources: &'a [MovingQuadrupoleSourceWindow<M>],
        history: &'a PreparedQuadrupoleHistory,
        microphone: FirstOrder,
        config: MovingPressureConfig,
        budget: MovingQuadrupoleFieldBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        let mut work = Work::new(budget.quadrature)?;
        if history.source_count() != sources.len() {
            return Err(sizing(
                "moving stress fields require the original record source count",
            ));
        }
        let [first_dd, last_dd] = history.source_time_bounds_s();
        let (first, last) = (first_dd.to_f64(), last_dd.to_f64());
        for (index, source) in sources.iter().enumerate() {
            check_cancel(cancelled)?;
            let initial = history
                .initial_source(index)
                .map_err(RenderError::Acoustic)?;
            if initial.stress_moment_n_m != [0.0; 6]
                || initial.stress_rate_n_m_s != [0.0; 6]
                || initial.stress_integral_n_m_s != [0.0; 6]
            {
                return Err(sizing(
                    "zero stress prehistory requires zero initial moment, rate and integral",
                ));
            }
            if source.source_time_interval_s[0] < first
                || source.source_time_interval_s[1] > last
                || source.source_bandwidth_hz * history.clock().timestep_s() >= 0.5
            {
                return Err(sizing(
                    "moving stress windows or bandwidth exceed the original record",
                ));
            }
            let position = source
                .motion
                .kinematics_at(first)
                .map_err(RenderError::Acoustic)?
                .position_m;
            if position != history.source_positions_m()[index] {
                return Err(sizing(
                    "moving stress source must start at its original record position",
                ));
            }
        }
        if config
            .max_block
            .checked_mul(size_of::<QuadrupoleFieldSample>())
            .is_none_or(|bytes| bytes > budget.pressure.maximum_scratch_bytes)
        {
            return Err(sizing(
                "moving stress field scratch exceeds its memory budget",
            ));
        }
        // Reuse the complete original pressure admission, including all image
        // roots, source windows and continuous Doppler/clearance bounds. Its
        // temporary scalar scratch is dropped before field scratch is allocated.
        let (path_count, maximum_retarded_iterations) = {
            let pressure = MovingQuadrupoleRenderer::new(
                observer,
                sources,
                |i, t| recorded_source(history, i, t),
                config,
                budget.pressure,
                cancelled,
            )?;
            (
                pressure.source_path_count(),
                pressure.maximum_retarded_iterations(),
            )
        };
        let medium = sources[0]
            .motion
            .observer_at(&observer, config.first_receiver_time_s)?
            .medium();
        let path_bytes = path_count
            .checked_mul(size_of::<Path>())
            .ok_or_else(|| sizing("moving stress path payload overflows"))?;
        let maximum_prefixes = budget
            .maximum_history_bytes
            .checked_sub(path_bytes)
            .ok_or_else(|| sizing("moving stress paths exceed the history payload budget"))?
            / size_of::<Prefix>();
        let mut paths = Vec::new();
        paths
            .try_reserve_exact(path_count)
            .map_err(|_| sizing("cannot reserve moving stress paths"))?;
        let last_receiver = super::clock_time(config, config.samples - 1)?;
        let mut prefix_count = 0;
        for (index, source) in sources.iter().enumerate() {
            for image in 0..source.motion.path_count() {
                check_cancel(cancelled)?;
                let final_path = source.motion.prepare_quadrupole_path_at(
                    image,
                    &observer,
                    last_receiver,
                    config.maximum_time_residual_s,
                    cancelled,
                )?;
                let receiver = source.motion.observer_at(&observer, last_receiver)?;
                let denominator = (1.0 - receiver.limits().maximum_mach).next_down();
                let padding =
                    ((2.0 * config.maximum_time_residual_s).next_up() / denominator).next_up();
                let covered_end = (final_path.source_time_s() + padding).next_up().min(last);
                let geometry = source.motion.force_history_path_at(
                    &observer,
                    image,
                    config.first_receiver_time_s,
                    [first, covered_end],
                    cancelled,
                )?;
                if !geometry.source_end_s.is_finite() || geometry.source_end_s < first {
                    return Err(sizing("moving stress history duration is unrepresentable"));
                }
                let start = prefix_count;
                prefix_count = charge(prefix_count, maximum_prefixes)?;
                let (mut time, mut knot) = (first, 0);
                while let Some((next, next_knot)) =
                    next_breakpoint(history, source.motion, time, knot, geometry.source_end_s)?
                {
                    check_cancel(cancelled)?;
                    prefix_count = charge(prefix_count, maximum_prefixes)?;
                    (time, knot) = (next, next_knot);
                }
                paths.push(Path {
                    source: index,
                    image,
                    geometry,
                    start,
                    count: prefix_count - start,
                });
            }
        }
        let mut prefixes = Vec::new();
        prefixes
            .try_reserve_exact(prefix_count)
            .map_err(|_| sizing("cannot reserve moving stress history prefixes"))?;
        for path in &paths {
            let source = &sources[path.source];
            let duration = path.geometry.source_end_s - first;
            let mut prefix = Prefix {
                time: first,
                source_knot: 0,
                integral: [Dd::ZERO; 3],
            };
            prefixes.push(prefix);
            while let Some((next, knot)) = next_breakpoint(
                history,
                source.motion,
                prefix.time,
                prefix.source_knot,
                path.geometry.source_end_s,
            )? {
                check_cancel(cancelled)?;
                let increment = if M::MOVING_HISTORY_KERNEL {
                    [Dd::ZERO; 3]
                } else {
                    integrate_cell(
                        history,
                        path.source,
                        source.motion,
                        path.geometry,
                        [prefix.time, next],
                        duration,
                        budget.quadrature,
                        &mut work,
                        cancelled,
                    )?
                };
                prefix = Prefix {
                    time: next,
                    source_knot: knot,
                    integral: add(prefix.integral, increment)?,
                };
                prefixes.push(prefix);
            }
        }
        let mut scratch = Vec::new();
        scratch
            .try_reserve_exact(config.max_block)
            .map_err(|_| sizing("cannot reserve moving stress field scratch"))?;
        scratch.resize(config.max_block, QuadrupoleFieldSample::default());
        check_cancel(cancelled)?;
        Ok(Self {
            observer,
            sources,
            history,
            config,
            microphone,
            medium,
            paths,
            prefixes,
            scratch,
            quadrature: budget.quadrature,
            work,
            completed: 0,
            maximum_retarded_iterations,
        })
    }

    pub const fn observer(&self) -> &M::Observer {
        &self.observer
    }
    pub const fn config(&self) -> MovingPressureConfig {
        self.config
    }
    pub const fn microphone(&self) -> FirstOrder {
        self.microphone
    }

    /// Background fluid velocity from the original shared laboratory frame.
    pub fn mean_flow_velocity_m_s(&self) -> [f64; 3] {
        self.sources[0].motion.mean_flow_velocity_m_s()
    }

    /// Energy and laboratory transport of one coherent stress-field sample.
    /// Retains the original material/frame and uses physical p/u independently
    /// of microphone projection. This observation changes no history or cursor.
    pub fn energy_transport(
        &self,
        field: &QuadrupoleFieldSample,
    ) -> Result<fs_aeroac::curle3d::energy::AcousticEnergyTransport, RenderError> {
        fs_aeroac::curle3d::energy::instantaneous(
            self.medium,
            self.mean_flow_velocity_m_s(),
            field.pressure_pa,
            field.particle_velocity_m_s,
        )
        .map_err(RenderError::Acoustic)
    }
    pub fn source_count(&self) -> usize {
        self.sources.len()
    }
    pub fn source_path_count(&self) -> usize {
        self.paths.len()
    }
    pub const fn samples_rendered(&self) -> u64 {
        self.completed
    }
    pub const fn max_block_len(&self) -> usize {
        self.config.max_block
    }
    pub const fn maximum_retarded_iterations(&self) -> u64 {
        self.maximum_retarded_iterations
    }
    pub const fn quadrature_evaluations(&self) -> u64 {
        self.work.used
    }
    pub fn remaining_samples(&self) -> u64 {
        self.config.samples - self.completed
    }
    pub fn receiver_time_at(&self, index: u64) -> Result<f64, RenderError> {
        if index >= self.config.samples {
            return Err(sizing(
                "moving stress sample exceeds the finite receiver window",
            ));
        }
        super::clock_time(self.config, index)
    }
    pub fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        if rate != self.config.sample_rate_hz {
            return Err(sizing("moving stress PCM and sink clocks differ"));
        }
        Ok(())
    }
    pub fn validate_sample_count(&self, count: u64) -> Result<(), RenderError> {
        if self
            .completed
            .checked_add(count)
            .is_none_or(|end| end > self.config.samples)
        {
            return Err(sizing(
                "moving stress PCM request exceeds the finite receiver window",
            ));
        }
        Ok(())
    }
    pub fn field_block(&mut self, output: &mut [QuadrupoleFieldSample]) -> Result<(), RenderError> {
        self.field_block_controlled(output, &mut || false)
    }
    pub fn field_block_controlled(
        &mut self,
        output: &mut [QuadrupoleFieldSample],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        self.prepare_block(output.len(), cancelled)?;
        output.copy_from_slice(&self.scratch[..output.len()]);
        self.completed += output.len() as u64;
        Ok(())
    }
    pub fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.block_controlled(output, &mut || false)
    }
    pub fn block_controlled(
        &mut self,
        output: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        self.prepare_block(output.len(), cancelled)?;
        for (target, field) in output.iter_mut().zip(&self.scratch) {
            *target = field.microphone_pressure_pa;
        }
        self.completed += output.len() as u64;
        Ok(())
    }

    fn prepare_block(
        &mut self,
        count: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        check_cancel(cancelled)?;
        if count == 0 {
            return Err(RenderError::EmptyBlock);
        }
        if count > self.config.max_block {
            return Err(sizing("moving stress field callback exceeds its capacity"));
        }
        self.validate_sample_count(count as u64)?;
        let first = self.history.clock().epoch_s();
        for i in 0..count {
            let time = self.receiver_time_at(self.completed + i as u64)?;
            let (mut pressure, mut velocity) = (Dd::ZERO, [Dd::ZERO; 3]);
            for path in &self.paths {
                check_cancel(cancelled)?;
                let source = &self.sources[path.source];
                let prepared = source.motion.prepare_quadrupole_field_path_at(
                    path.image,
                    &self.observer,
                    time,
                    self.config.maximum_time_residual_s,
                    cancelled,
                )?;
                let emission = prepared.source_time_s();
                let bank = &self.prefixes[path.start..path.start + path.count];
                let upper = bank.partition_point(|prefix| prefix.time <= emission);
                if upper == 0
                    || emission > path.geometry.source_end_s
                    || emission < source.source_time_interval_s[0]
                    || emission > source.source_time_interval_s[1]
                {
                    return Err(sizing(
                        "moving stress query exceeds its original source or past geometry",
                    ));
                }
                let geometry = if M::MOVING_HISTORY_KERNEL {
                    source.motion.force_history_path_at(
                        &self.observer,
                        path.image,
                        time,
                        [first, path.geometry.source_end_s],
                        cancelled,
                    )?
                } else {
                    path.geometry
                };
                let duration = geometry.source_end_s - first;
                let prefix = bank[upper - 1];
                let mut history_integral = prefix.integral;
                if M::MOVING_HISTORY_KERNEL {
                    for cell in bank[..upper].windows(2) {
                        let increment = integrate_cell(
                            self.history,
                            path.source,
                            source.motion,
                            geometry,
                            [cell[0].time, cell[1].time],
                            duration,
                            self.quadrature,
                            &mut self.work,
                            cancelled,
                        )?;
                        history_integral = add(history_integral, increment)?;
                    }
                }
                let increment = integrate_cell(
                    self.history,
                    path.source,
                    source.motion,
                    geometry,
                    [prefix.time, emission],
                    duration,
                    self.quadrature,
                    &mut self.work,
                    cancelled,
                )?;
                history_integral = add(history_integral, increment)?;
                let strength = recorded_source(self.history, path.source, emission)
                    .map_err(RenderError::Acoustic)?;
                // The integral uses the untransformed physical tensor at the
                // virtual receiver. This field owner applies image parity and
                // the polar-vector transform to local AND history terms once.
                let (p, u) = prepared
                    .evaluate_dd(strength, history_integral)
                    .map_err(RenderError::Acoustic)?;
                pressure = pressure + p;
                for axis in 0..3 {
                    velocity[axis] = velocity[axis] + u[axis];
                }
            }
            self.scratch[i] = QuadrupoleFieldSample::observe(
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

impl<M: QuadrupoleMotion + ForceMotion> PressureRenderer for MovingQuadrupoleFieldRenderer<'_, M> {
    fn samples_rendered(&self) -> u64 {
        Self::samples_rendered(self)
    }
    fn max_block_len(&self) -> usize {
        Self::max_block_len(self)
    }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        Self::validate_sample_rate(self, rate)
    }
    fn validate_sample_count(&self, count: u64) -> Result<(), RenderError> {
        Self::validate_sample_count(self, count)
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        Self::block(self, output)
    }
}

fn recorded_source(
    history: &PreparedQuadrupoleHistory,
    source: usize,
    time: f64,
) -> Result<TransientQuadrupoleSource, AeroacError> {
    let [first, last] = history.source_time_bounds_s();
    // Preserve an original endpoint whose DD timestamp rounds to this f64
    // motion clock. This is exact endpoint identification, not extrapolation.
    let query = if time == first.to_f64() {
        first
    } else if time == last.to_f64() {
        last
    } else {
        Dd::from_f64(time)
    };
    history.source_at_clock_dd(source, query)
}

fn charge(count: usize, maximum: usize) -> Result<usize, RenderError> {
    count
        .checked_add(1)
        .filter(|next| *next <= maximum)
        .ok_or_else(|| sizing("moving stress prefixes exceed the history payload budget"))
}

fn next_breakpoint<M: ForceMotion>(
    history: &PreparedQuadrupoleHistory,
    motion: M,
    time: f64,
    source_knot: usize,
    end: f64,
) -> Result<Option<(f64, usize)>, RenderError> {
    if time >= end {
        return Ok(None);
    }
    let next_source = if source_knot + 1 < history.sample_count() {
        history
            .clock()
            .time_at(source_knot + 1)
            .map_err(RenderError::Acoustic)?
    } else {
        f64::INFINITY
    };
    let next_geometry = motion
        .next_force_history_knot_after(time)?
        .unwrap_or(f64::INFINITY);
    let next = next_source.min(next_geometry);
    if next > end {
        return Ok(None);
    }
    if !next.is_finite() || next <= time {
        return Err(sizing("moving stress prefix clock lost its original order"));
    }
    Ok(Some((next, source_knot + usize::from(next == next_source))))
}

fn add(first: [Dd; 3], second: [Dd; 3]) -> Result<[Dd; 3], RenderError> {
    let result: [Dd; 3] = core::array::from_fn(|i| first[i] + second[i]);
    if result.iter().any(|value| !value.to_f64().is_finite()) {
        return Err(sizing("moving stress history integral is unrepresentable"));
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn integrate_cell<M: ForceMotion>(
    history: &PreparedQuadrupoleHistory,
    source_index: usize,
    motion: M,
    geometry: ForceHistoryPath,
    interval: [f64; 2],
    duration: f64,
    quadrature: ForceHistoryQuadrature,
    work: &mut Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Dd; 3], RenderError> {
    if interval[0] == interval[1] {
        check_cancel(cancelled)?;
        return Ok([Dd::ZERO; 3]);
    }
    let tolerance = quadrature.absolute_tolerance_n_s_m3 * ((interval[1] - interval[0]) / duration);
    integral::integrate(
        interval[0],
        interval[1],
        tolerance,
        geometry.maximum_source_speed_m_s,
        geometry.minimum_history_distance_m,
        quadrature,
        work,
        |time| {
            let strength =
                recorded_source(history, source_index, time).map_err(RenderError::Acoustic)?;
            let position = motion
                .kinematics_at(time)
                .map_err(RenderError::Acoustic)?
                .position_m;
            let receiver = geometry.observer.observer_position_m();
            let range = det::hypot(
                det::hypot(receiver[0] - position[0], receiver[1] - position[1]),
                receiver[2] - position[2],
            );
            let value = geometry
                .observer
                .quadrupole_history_integrand(position, strength.stress_moment_n_m)
                .map_err(RenderError::Acoustic)?;
            Ok((value, range))
        },
        cancelled,
    )
}
