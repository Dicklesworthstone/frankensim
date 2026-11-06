//! Causal volume-stress histories through the existing microphone/PCM owner.
//!
//! The original transient quadrupole observer owns propagation, including each
//! reflected path's delay and tensor transformation. This adapter samples the
//! original physical source history at those emission times, projects the full
//! coherent pressure/velocity field onto the native ideal microphone, and
//! implements the same `PressureRenderer` used by decimation and WAV streaming.
//! No mechanics advances, source interpolation, delay compensation, artificial
//! tail, normalizing gain or second audio encoder is introduced here.

use super::PressureRenderer;
use crate::render::RenderError;
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::quadrupole::transient::{
    PreparedTransientQuadrupoleObserver, TransientQuadrupoleSource,
};
use fs_aeroac::observation::HarmonicSampleClock;
use fs_bem::near_field::FirstOrder;
use fs_math::dd::Dd;

/// Doppler-correct pressure from uniform and accelerating volume-stress histories.
pub mod moving;

/// Physical pressure, particle motion, intensity and ideal microphone output.
pub type QuadrupoleFieldSample = super::moving::field::CompactFieldSample;

/// Declared support of one immutable physical source's history.
///
/// The source reader must return its consistent Q, Qdot, Qddot and integral
/// throughout this closed interval. Its bandwidth includes the actual stress
/// history; numerical interpolation is not an antialias filter. An initial
/// stress integral belongs to that history and never resets at audio blocks.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadrupoleSourceWindow {
    /// First and last available emission times [s] on the original source clock.
    pub source_time_interval_s: [Dd; 2],
    /// Declared source bandwidth [Hz]; zero admits a constant source.
    pub source_bandwidth_hz: f64,
}

/// One finite receiver clock, before optional causal pressure decimation.
#[derive(Clone, Copy, Debug)]
pub struct QuadrupolePressureConfig {
    /// Native uniform clock with a double-double epoch; no harmonic source
    /// assumption is added by reusing its exact integer-index timestamps.
    pub clock: HarmonicSampleClock,
    /// Total receiver samples in the finite observation window.
    pub samples: u64,
    /// Largest admitted callback, in receiver samples.
    pub max_block: usize,
    /// Must cover every stationary source and be strictly below Nyquist.
    pub receiver_bandwidth_hz: f64,
}

/// Bounds on this observation, separate from the original history/preparation.
#[derive(Clone, Copy, Debug)]
pub struct QuadrupolePressureBudget {
    /// Physical sources, before the observer expands reflected image paths.
    pub maximum_sources: usize,
    /// Total samples on the original receiver clock, before decimation.
    pub maximum_output_samples: u64,
    /// All physical source/image paths times the complete output window.
    pub maximum_source_path_samples: u64,
    /// Owned field scratch only; borrowed observers and source histories retain
    /// their original preparation budgets. Allocator overhead is excluded.
    pub maximum_scratch_bytes: usize,
}

/// Reflected transient volume stresses on the native pressure streaming path.
///
/// `source_at(index, emission_time)` is a read-only query of the original
/// physical history, in the observer's original source order and frame. The
/// SAME source can be queried repeatedly at different times for its images.
/// The full receiver window is admitted against every path before allocation.
/// Preparation owns one bounded scratch buffer. Callbacks allocate nothing,
/// retain their global integer clock, and publish output/cursor atomically.
/// A source reader with hidden side effects violates this observation contract.
/// Emission times retain DD precision: subtract the source's epoch or evaluate
/// its phase before rounding. A field sample's f64 timestamp is for inspection;
/// `receiver_time_at` returns the original full-precision clock.
///
/// The observer's ideal stationary homogeneous lossless propagation scope is
/// retained. Neither declared bandwidth nor successful PCM rendering validates
/// source compactness, a turbulence model or interpolation convergence.
pub struct TransientQuadrupoleRenderer<'a, F> {
    observer: &'a PreparedTransientQuadrupoleObserver,
    source_windows: &'a [QuadrupoleSourceWindow],
    source_at: F,
    microphone: FirstOrder,
    config: QuadrupolePressureConfig,
    scratch: Vec<QuadrupoleFieldSample>,
    completed: u64,
}

impl<'a, F> TransientQuadrupoleRenderer<'a, F>
where
    F: Fn(usize, Dd) -> Result<TransientQuadrupoleSource, AeroacError>,
{
    pub fn new(
        observer: &'a PreparedTransientQuadrupoleObserver,
        source_windows: &'a [QuadrupoleSourceWindow],
        source_at: F,
        microphone: FirstOrder,
        config: QuadrupolePressureConfig,
        budget: QuadrupolePressureBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        let rate = f64::from(config.clock.sample_rate_hz());
        if observer.source_count() == 0
            || observer.source_count() != source_windows.len()
            || source_windows.len() > budget.maximum_sources
            || config.samples == 0
            || config.samples > budget.maximum_output_samples
            || config.max_block == 0
            || !config.receiver_bandwidth_hz.is_finite()
            || config.receiver_bandwidth_hz < 0.0
            || config.receiver_bandwidth_hz >= 0.5 * rate
        {
            return Err(sizing(
                "quadrupole PCM requires finite clocks, bandwidth and capacities",
            ));
        }
        if u64::try_from(observer.path_count())
            .ok()
            .and_then(|paths| paths.checked_mul(config.samples))
            .is_none_or(|work| work > budget.maximum_source_path_samples)
        {
            return Err(sizing(
                "quadrupole PCM exceeds its complete source/image work budget",
            ));
        }
        if config
            .max_block
            .checked_mul(size_of::<QuadrupoleFieldSample>())
            .is_none_or(|bytes| bytes > budget.maximum_scratch_bytes)
        {
            return Err(sizing("quadrupole PCM scratch exceeds its memory budget"));
        }
        let first = clock_time(config, 0)?;
        let last = clock_time(config, config.samples - 1)?;
        let step = Dd::ONE / Dd::from_f64(rate);
        if !first.lt(checked_time(first + step)?) || !last.lt(checked_time(last + step)?) {
            return Err(sizing(
                "quadrupole receiver sample spacing is unrepresentable",
            ));
        }
        for (index, source) in source_windows.iter().enumerate() {
            check_cancel(cancelled)?;
            let [start, end] = source.source_time_interval_s;
            let start = checked_time(start)?;
            let end = checked_time(end)?;
            if !start.lt(end)
                || !source.source_bandwidth_hz.is_finite()
                || source.source_bandwidth_hz < 0.0
                || source.source_bandwidth_hz > config.receiver_bandwidth_hz
            {
                return Err(sizing(
                    "quadrupole histories require finite support and resolved bandwidth",
                ));
            }
            let [shortest, longest] = observer
                .source_delay_bounds_s(index)
                .ok_or_else(|| sizing("quadrupole source has no admitted propagation paths"))?;
            let earliest = retarded_time(first, longest)?;
            let latest = retarded_time(last, shortest)?;
            // The shortest delay must remain resolved even at the first end
            // of a window whose largest absolute timestamp occurs there.
            retarded_time(first, shortest)?;
            if earliest.lt(start) || end.lt(latest) {
                return Err(sizing(
                    "quadrupole receiver window lacks complete direct/reflected history",
                ));
            }
        }
        let mut scratch = Vec::new();
        scratch
            .try_reserve_exact(config.max_block)
            .map_err(|_| sizing("cannot reserve quadrupole PCM scratch"))?;
        scratch.resize(config.max_block, QuadrupoleFieldSample::default());
        check_cancel(cancelled)?;
        Ok(Self {
            observer,
            source_windows,
            source_at,
            microphone,
            config,
            scratch,
            completed: 0,
        })
    }

    pub const fn observer(&self) -> &PreparedTransientQuadrupoleObserver {
        self.observer
    }
    pub const fn source_windows(&self) -> &[QuadrupoleSourceWindow] {
        self.source_windows
    }
    pub const fn microphone(&self) -> FirstOrder {
        self.microphone
    }
    pub const fn config(&self) -> QuadrupolePressureConfig {
        self.config
    }
    pub const fn samples_rendered(&self) -> u64 {
        self.completed
    }
    pub const fn max_block_len(&self) -> usize {
        self.config.max_block
    }

    /// Global receiver timestamp; block partitioning never changes arithmetic.
    pub fn receiver_time_at(&self, index: u64) -> Result<Dd, RenderError> {
        if index >= self.config.samples {
            return Err(sizing(
                "quadrupole receiver timestamp exceeds its finite window",
            ));
        }
        clock_time(self.config, index)
    }

    pub fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        if rate != self.config.clock.sample_rate_hz() {
            return Err(sizing("quadrupole PCM and sink sample clocks differ"));
        }
        Ok(())
    }

    pub fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self
            .completed
            .checked_add(samples)
            .is_none_or(|end| end > self.config.samples)
        {
            return Err(sizing(
                "quadrupole PCM request exceeds the finite history window",
            ));
        }
        Ok(())
    }

    /// Publish full fields; they retain the original observation rate.
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

    /// Publish Pa-equivalent microphone pressure to the shared PCM/decimator.
    pub fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.block_controlled(output, &mut || false)
    }

    pub fn block_controlled(
        &mut self,
        output: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        self.prepare_block(output.len(), cancelled)?;
        for (out, field) in output.iter_mut().zip(&self.scratch) {
            *out = field.microphone_pressure_pa;
        }
        self.completed += output.len() as u64;
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
        if count > self.config.max_block {
            return Err(sizing("quadrupole PCM block exceeds its callback capacity"));
        }
        self.validate_sample_count(
            u64::try_from(count).map_err(|_| sizing("quadrupole block count overflow"))?,
        )?;
        for i in 0..count {
            check_cancel(cancelled)?;
            let time = self.receiver_time_at(self.completed + i as u64)?;
            let field = self
                .observer
                .evaluate_at_clock_dd_cancellable(
                    time,
                    |source, emission| (self.source_at)(source, emission),
                    &mut *cancelled,
                )
                .map_err(RenderError::Acoustic)?;
            self.scratch[i] = QuadrupoleFieldSample::observe(
                time.to_f64(),
                Dd::from_f64(field.pressure_pa),
                field.particle_velocity_m_s.map(Dd::from_f64),
                self.microphone,
                self.observer.acoustic(),
            )?;
        }
        check_cancel(cancelled)
    }
}

impl<F> PressureRenderer for TransientQuadrupoleRenderer<'_, F>
where
    F: Fn(usize, Dd) -> Result<TransientQuadrupoleSource, AeroacError>,
{
    fn samples_rendered(&self) -> u64 {
        Self::samples_rendered(self)
    }
    fn max_block_len(&self) -> usize {
        Self::max_block_len(self)
    }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        Self::validate_sample_rate(self, rate)
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        Self::validate_sample_count(self, samples)
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        Self::block(self, output)
    }
}

fn clock_time(config: QuadrupolePressureConfig, index: u64) -> Result<Dd, RenderError> {
    checked_time(config.clock.time_at(index).map_err(RenderError::Acoustic)?)
}

fn retarded_time(receiver: Dd, delay: f64) -> Result<Dd, RenderError> {
    let time = checked_time(receiver - Dd::from_f64(delay))?;
    if !time.lt(receiver) {
        return Err(sizing("quadrupole retarded clock is unrepresentable"));
    }
    Ok(time)
}

fn checked_time(time: Dd) -> Result<Dd, RenderError> {
    if !time.hi.is_finite() || !time.lo.is_finite() {
        return Err(sizing("quadrupole source/receiver clock is nonfinite"));
    }
    let normalized = Dd::from_pair(time.hi, time.lo);
    if !normalized.hi.is_finite() || !normalized.lo.is_finite() {
        return Err(sizing("quadrupole source/receiver clock overflows"));
    }
    Ok(normalized)
}

fn sizing(what: &'static str) -> RenderError {
    RenderError::Sizing { what }
}

fn check_cancel(cancelled: &mut impl FnMut() -> bool) -> Result<(), RenderError> {
    if cancelled() {
        return Err(RenderError::Acoustic(AeroacError::Cancelled {
            what: "transient quadrupole PCM",
        }));
    }
    Ok(())
}
