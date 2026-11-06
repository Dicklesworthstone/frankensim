//! Physical moving stress histories on the shared pressure/PCM streaming path.
//!
//! Each callback queries Q, Qdot and Qddot on its original emission clock and
//! evaluates the complete retarded quadrupole pressure. The accepted tensor is
//! in fixed axes of the homogeneous medium rest frame, including its trace;
//! translation neither rotates nor renormalizes it. No history interpolation,
//! source evolution or artificial prehistory occurs here. Received bandwidth
//! must cover Doppler shift AND geometric amplitude modulation. Render at an
//! adequate rate and use the existing `DecimatedRenderer` for downsampling.
//!
//! Scope: uniform or original admitted quintic source translations, a fixed
//! omnidirectional receiver, and homogeneous quiescent lossless propagation.
//! Analytic trajectory jerk supplies the accelerating radiation coefficients.
//! C2 knots can have jerk jumps: the original clock-rate bound does not certify
//! the bandwidth of that geometric modulation. This pressure-only
//! adapter does not infer particle velocity, directivity from p/Z, reflected
//! moving paths, compactness, or the validity of a turbulence source model.

use super::{QuadrupolePressureBudget, check_cancel, sizing};
use crate::pcm_wav::observation::PressureRenderer;
use crate::pcm_wav::observation::moving::{MovingPressureConfig, clock_time};
use crate::render::RenderError;
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::{MovingCompactObserver, UniformSourceMotion};
use fs_aeroac::curle3d::quadrupole::transient::TransientQuadrupoleSource;
use fs_math::dd::Dd;

mod motion;
pub use motion::QuadrupoleMotion;

/// One original physical stress history and its independent admitted motion.
#[derive(Clone, Copy, Debug)]
pub struct MovingQuadrupoleSourceWindow<M = UniformSourceMotion> {
    /// Translation in the same medium rest frame as the receiver and tensor.
    pub motion: M,
    /// Closed interval [s] available from the original source reader.
    pub source_time_interval_s: [f64; 2],
    /// Physical history bandwidth [Hz], before Doppler and amplitude modulation.
    pub source_bandwidth_hz: f64,
}

/// Bounded, coherent Doppler pressure feeding the existing decimator and WAV.
///
/// `source_at(index, emission_time)` must be a read-only query of the same
/// immutable physical source throughout its declared interval, returning
/// consistent Q/Qdot/Qddot in medium-frame axes. Its causal stress integral is
/// unused by pressure. Producer interpolation and convergence remain explicit
/// responsibilities. Times retain the existing moving operator's f64 clock
/// domain; clocks too large to resolve adjacent samples are refused.
///
/// Admission covers the entire observation window, including continuous
/// closest approach and the maximum retarded clock rate. Scratch is allocated
/// once. Each block uses global integer indices and fixed-order DD source sums,
/// allocates nothing, and publishes output/cursor only after the whole block
/// succeeds. Failed reads or cancellation leave both unchanged.
pub struct MovingQuadrupoleRenderer<'a, F, M: QuadrupoleMotion = UniformSourceMotion> {
    observer: MovingCompactObserver,
    sources: &'a [MovingQuadrupoleSourceWindow<M>],
    source_at: F,
    config: MovingPressureConfig,
    scratch: Vec<f64>,
    completed: u64,
    maximum_retarded_iterations: u64,
}

impl<'a, F, M: QuadrupoleMotion> MovingQuadrupoleRenderer<'a, F, M>
where
    F: Fn(usize, f64) -> Result<TransientQuadrupoleSource, AeroacError>,
{
    /// Admit all clocks/history/geometry/bandwidth and finite lifetime budgets.
    pub fn new(
        observer: MovingCompactObserver,
        sources: &'a [MovingQuadrupoleSourceWindow<M>],
        source_at: F,
        config: MovingPressureConfig,
        budget: QuadrupolePressureBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RenderError> {
        check_cancel(cancelled)?;
        let rate = f64::from(config.sample_rate_hz);
        if sources.is_empty()
            || sources.len() > budget.maximum_sources
            || config.sample_rate_hz == 0
            || config.samples == 0
            || config.samples > budget.maximum_output_samples
            || config.samples > (1_u64 << 53)
            || config.max_block == 0
            || !config.first_receiver_time_s.is_finite()
            || !config.receiver_bandwidth_hz.is_finite()
            || config.receiver_bandwidth_hz < 0.0
            || config.receiver_bandwidth_hz >= 0.5 * rate
            || !config.maximum_time_residual_s.is_finite()
            || !(0.0..=0.25 / rate).contains(&config.maximum_time_residual_s)
        {
            return Err(sizing(
                "moving quadrupole requires finite clocks, bandwidth and capacities",
            ));
        }
        let work = (sources.len() as u64)
            .checked_mul(config.samples)
            .ok_or_else(|| sizing("moving quadrupole source/sample work overflows"))?;
        let bytes = config
            .max_block
            .checked_mul(core::mem::size_of::<f64>())
            .ok_or_else(|| sizing("moving quadrupole scratch sizing overflows"))?;
        if work > budget.maximum_source_path_samples || bytes > budget.maximum_scratch_bytes {
            return Err(sizing("moving quadrupole exceeds work or scratch budget"));
        }
        let first_time = clock_time(config, 0)?;
        let last_time = clock_time(config, config.samples - 1)?;
        if config.samples > 1
            && (clock_time(config, 1)? <= first_time
                || clock_time(config, config.samples - 2)? >= last_time)
        {
            return Err(sizing(
                "moving quadrupole receiver clock cannot resolve adjacent samples",
            ));
        }
        let mut iterations_per_sample = 0_u64;
        for source in sources {
            check_cancel(cancelled)?;
            source.motion.admit_observer(&observer)?;
            source
                .motion
                .admit_clock(first_time, last_time, config.maximum_time_residual_s)?;
            iterations_per_sample = iterations_per_sample
                .checked_add(
                    u64::try_from(source.motion.maximum_root_iterations())
                        .map_err(|_| sizing("moving quadrupole root limit exceeds u64"))?,
                )
                .ok_or_else(|| sizing("moving quadrupole root work overflows"))?;
            let [start, end] = source.source_time_interval_s;
            if !start.is_finite()
                || !end.is_finite()
                || start > end
                || !source.source_bandwidth_hz.is_finite()
                || source.source_bandwidth_hz < 0.0
            {
                return Err(sizing(
                    "moving quadrupole needs finite ordered source support and bandwidth",
                ));
            }
            let first = source.motion.prepare_quadrupole_at(
                &observer,
                first_time,
                config.maximum_time_residual_s,
                cancelled,
            )?;
            let last = source.motion.prepare_quadrupole_at(
                &observer,
                last_time,
                config.maximum_time_residual_s,
                cancelled,
            )?;
            if first.source_time_s() < start
                || last.source_time_s() > end
                || first.source_time_s() > last.source_time_s()
            {
                return Err(sizing(
                    "moving quadrupole retarded clock exceeds declared source support",
                ));
            }
            // Each original motion owner admits continuous clearance and
            // bounds d tau/d t over this whole finite receiver window.
            let maximum_rate = source.motion.admit_window(
                &observer,
                first.propagation(),
                last.propagation(),
                config.maximum_time_residual_s,
                cancelled,
            )?;
            let shifted_bandwidth = source.source_bandwidth_hz * maximum_rate;
            if !shifted_bandwidth.is_finite() || shifted_bandwidth > config.receiver_bandwidth_hz {
                return Err(sizing(
                    "moving quadrupole received bandwidth does not cover Doppler shift",
                ));
            }
        }
        let maximum_retarded_iterations = iterations_per_sample
            .checked_mul(config.samples)
            .ok_or_else(|| sizing("moving quadrupole lifetime root work overflows"))?;
        check_cancel(cancelled)?;
        let mut scratch = Vec::new();
        scratch
            .try_reserve_exact(config.max_block)
            .map_err(|_| sizing("moving quadrupole scratch allocation refused"))?;
        scratch.resize(config.max_block, 0.0);
        Ok(Self {
            observer,
            sources,
            source_at,
            config,
            scratch,
            completed: 0,
            maximum_retarded_iterations,
        })
    }

    pub const fn samples_rendered(&self) -> u64 {
        self.completed
    }

    pub const fn max_block_len(&self) -> usize {
        self.config.max_block
    }

    /// Complete-window cap on interior root refinements, excluding admission
    /// and the constant-work endpoint/final evaluations. Uniform roots use 0.
    pub const fn maximum_retarded_iterations(&self) -> u64 {
        self.maximum_retarded_iterations
    }

    pub fn receiver_time_at(&self, index: u64) -> Result<f64, RenderError> {
        if index >= self.config.samples {
            return Err(sizing(
                "moving quadrupole receiver timestamp exceeds its finite window",
            ));
        }
        clock_time(self.config, index)
    }

    pub fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        if rate != self.config.sample_rate_hz {
            return Err(sizing("moving quadrupole PCM and sink clocks differ"));
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
                "moving quadrupole PCM request exceeds its finite window",
            ));
        }
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
        check_cancel(cancelled)?;
        if output.is_empty() || output.len() > self.config.max_block {
            return Err(sizing(
                "moving quadrupole callback is empty or exceeds its bound",
            ));
        }
        self.validate_sample_count(output.len() as u64)?;
        for i in 0..output.len() {
            let time = clock_time(self.config, self.completed + i as u64)?;
            let mut sum = Dd::ZERO;
            for (index, source) in self.sources.iter().enumerate() {
                check_cancel(cancelled)?;
                let prepared = source.motion.prepare_quadrupole_at(
                    &self.observer,
                    time,
                    self.config.maximum_time_residual_s,
                    cancelled,
                )?;
                let source_time = prepared.source_time_s();
                // Keep source readers inside the original finite support even
                // if rounded retarded roots sit just outside an endpoint.
                if source_time < source.source_time_interval_s[0]
                    || source_time > source.source_time_interval_s[1]
                {
                    return Err(sizing(
                        "moving quadrupole source query exceeds declared support",
                    ));
                }
                let strength =
                    (self.source_at)(index, source_time).map_err(RenderError::Acoustic)?;
                let pressure = prepared
                    .pressure_pa(strength)
                    .map_err(RenderError::Acoustic)?;
                sum = sum + Dd::from_f64(pressure);
            }
            let pressure = sum.to_f64();
            if !pressure.is_finite() {
                return Err(sizing(
                    "moving quadrupole coherent pressure is unrepresentable",
                ));
            }
            self.scratch[i] = pressure;
        }
        check_cancel(cancelled)?;
        output.copy_from_slice(&self.scratch[..output.len()]);
        self.completed += output.len() as u64;
        Ok(())
    }
}

impl<F, M: QuadrupoleMotion> PressureRenderer for MovingQuadrupoleRenderer<'_, F, M>
where
    F: Fn(usize, f64) -> Result<TransientQuadrupoleSource, AeroacError>,
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
