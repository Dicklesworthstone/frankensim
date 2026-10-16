//! Pressure observation on a different clock from the mechanics.
//!
//! Mechanics and every event stay on the admitted source clock. Only pressure
//! passes through the existing causal integer decimator; forces, velocities,
//! constitutive laws and state budgets are never filtered or retuned. Delay is
//! retained explicitly. This does not prove convergence of the high-rate model
//! or suppress frequencies already aliased at its own mechanical sample rate.

use super::decimate::Decimator;
use crate::render::{RenderError, schedule::ScheduledRenderer};

/// Latency-aligned acoustic superposition of independently clocked mechanics.
pub mod ensemble;

/// A mono pressure producer on its OUTPUT sample clock.
///
/// Admission methods must not advance any state. A failed physical callback
/// can leave a prefix advanced; implementations must then refuse continuation.
/// The stream writer owns cancellation boundaries and the PCM conversion.
pub trait PressureRenderer {
    /// Complete output samples, not internal mechanics steps or filter stages.
    fn samples_rendered(&self) -> u64;
    /// Largest admitted output callback in samples.
    fn max_block_len(&self) -> usize;
    /// Check both the output rate and the producer's healthy state.
    fn validate_sample_rate(&self, sample_rate_hz: u32) -> Result<(), RenderError>;
    /// Check all clocks for a whole request, before its first callback.
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError>;
    /// Produce a complete nonempty pressure callback within the admitted size.
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError>;
}

// Type erasure is optional and occurs at composition, never per sample. This
// lets a nonlinear impact observer and a modal/reed renderer share an output
// without forcing either to adopt the other's physical integration scheme.
impl<T: PressureRenderer + ?Sized> PressureRenderer for Box<T> {
    fn samples_rendered(&self) -> u64 { (**self).samples_rendered() }
    fn max_block_len(&self) -> usize { (**self).max_block_len() }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        (**self).validate_sample_rate(rate)
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        (**self).validate_sample_count(samples)
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        (**self).block(output)
    }
}

impl PressureRenderer for ScheduledRenderer {
    fn samples_rendered(&self) -> u64 { ScheduledRenderer::samples_rendered(self) }
    fn max_block_len(&self) -> usize { self.context().max_block_len() }
    fn validate_sample_rate(&self, sample_rate_hz: u32) -> Result<(), RenderError> {
        ScheduledRenderer::validate_sample_rate(self, sample_rate_hz)
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        ScheduledRenderer::validate_sample_count(self, samples)
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        ScheduledRenderer::block(self, output)
    }
}

impl PressureRenderer for crate::render::schedule::force::ensemble::EnsembleRenderer {
    fn samples_rendered(&self) -> u64 { Self::samples_rendered(self) }
    fn max_block_len(&self) -> usize { self.renderer().context().max_block_len() }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.renderer().validate_sample_rate(rate)
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        self.renderer().validate_sample_count(samples)?;
        if samples > self.remaining_samples() {
            return Err(sizing("requested pressure window exceeds the finite ensemble horizon"));
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        Self::block(self, output)
    }
}

/// Immutable observation metadata; a profile is not a physical certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecimationInfo {
    /// Actual rate of the already-constructed physical voices [Hz].
    pub mechanics_sample_rate_hz: u32,
    /// Rate of emitted pressure samples [Hz].
    pub output_sample_rate_hz: u32,
    /// Number of mechanics samples per output sample, in 1..=16.
    pub ratio: usize,
    /// Causal group delay, expressed in output samples. It is NOT removed.
    pub delay_output_samples: f64,
    /// Last source-sample index of the first emitted group: ratio - 1.
    pub first_output_source_index: usize,
    /// Fixed coefficient-design identity used by the shared decimator.
    pub filter_profile: &'static str,
}

/// An admitted pressure producer observed through the shared causal decimator.
/// ScheduledRenderer remains the default type for existing callers; nonlinear
/// impact or other physical producers use the SAME filter and streaming path.
///
/// Starts at source sample zero with zero filter history. Constructing a new
/// filter over an already-running source is refused rather than discarding its
/// missing pressure history. Continue by retaining this entire object.
///
/// Only one mechanics group is buffered, independent of duration or callback
/// size. Output errors poison the wrapper because the mechanical group or an
/// earlier output prefix may already have advanced. No rollback of a failed
/// callback, hard-real-time guarantee, or allocation-free voice is inferred.
pub struct DecimatedRenderer<S = ScheduledRenderer> {
    source: S,
    filter: Decimator,
    input: Vec<f64>,
    info: DecimationInfo,
    max_block: usize,
    completed: u64,
    poisoned: bool,
}
impl<S: PressureRenderer> DecimatedRenderer<S> {
    /// Bind an integer output clock to existing mechanics WITHOUT rebuilding it.
    /// Source controls retain their exact source-sample indices, including
    /// assignments between output boundaries. Noninteger ratios and upsampling
    /// refuse. Ratio one delegates the original callback without extra arithmetic.
    pub fn new(source: S, mechanics_sample_rate_hz: u32,
        output_sample_rate_hz: u32, max_block: usize) -> Result<Self, RenderError>
    {
        source.validate_sample_rate(mechanics_sample_rate_hz)?;
        if output_sample_rate_hz == 0 || mechanics_sample_rate_hz % output_sample_rate_hz != 0 {
            return Err(sizing("observation requires a positive integer mechanics/output rate ratio"));
        }
        let ratio = (mechanics_sample_rate_hz / output_sample_rate_hz) as usize;
        if !(1..=16).contains(&ratio) || max_block == 0 || source.max_block_len() == 0 {
            return Err(sizing("observation requires a ratio in 1..=16 and positive callback capacities"));
        }
        if source.samples_rendered() != 0 {
            return Err(sizing("decimation must start at source sample zero with complete filter history"));
        }
        if ratio == 1 && max_block > source.max_block_len() {
            return Err(sizing("bypass output capacity must fit the source callback"));
        }
        let filter = Decimator::new(ratio, 1).map_err(|what| RenderError::Control { what })?;
        let info = DecimationInfo {
            mechanics_sample_rate_hz, output_sample_rate_hz, ratio,
            delay_output_samples: filter.delay_output_frames(), first_output_source_index: ratio - 1,
            filter_profile: if ratio == 1 { "identity" } else { "blackman-harris-integer-decimator-v1" },
        };
        Ok(Self { source, filter, input: vec![0.0; ratio], info, max_block,
            completed: 0, poisoned: false })
    }

    /// Complete observation policy, including uncompensated causal delay.
    #[must_use]
    pub const fn info(&self) -> DecimationInfo { self.info }
    /// Read source state without bypassing the retained observation history.
    #[must_use]
    pub const fn source(&self) -> &S { &self.source }
    /// Number of complete OUTPUT samples. Discard output from a failed callback.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 { self.completed }
    /// Admitted output callback capacity, independent of the source capacity.
    #[must_use]
    pub const fn max_block_len(&self) -> usize { self.max_block }
    /// No rate conversion is silently selected by a sink.
    pub fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        self.source.validate_sample_rate(self.info.mechanics_sample_rate_hz)?;
        if rate != self.info.output_sample_rate_hz {
            return Err(sizing("sink rate must match the declared decimated output clock"));
        }
        Ok(())
    }
    /// Admit the complete output request and its expanded mechanical work clock.
    pub fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        self.completed.checked_add(samples).ok_or_else(|| sizing("output sample clock overflow"))?;
        let mechanics = samples.checked_mul(self.info.ratio as u64)
            .ok_or_else(|| sizing("expanded mechanics sample count overflow"))?;
        PressureRenderer::validate_sample_count(&self.source, mechanics)
    }
    /// Convert an exact mechanics horizon without padding, truncation or flushing.
    /// A fractional last output interval refuses, even when no event occurs there.
    pub fn output_samples_for(&self, mechanics_samples: u64) -> Result<u64, RenderError> {
        let ratio = self.info.ratio as u64;
        if mechanics_samples % ratio != 0 {
            return Err(sizing("mechanics horizon must contain complete output intervals; no samples are dropped"));
        }
        let samples = mechanics_samples / ratio;
        self.validate_sample_count(samples)?;
        Ok(samples)
    }
    /// Render one complete output callback. Its internal mechanics groups are
    /// not cancellation boundaries; an enclosing gated stream drains the callback.
    /// No synthetic tail samples are produced and group delay is not compensated.
    pub fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.max_block { return Err(sizing("output block exceeds decimated callback capacity")); }
        let samples = u64::try_from(output.len()).map_err(|_| sizing("output length exceeds u64"))?;
        self.validate_sample_count(samples)?;
        if self.info.ratio == 1 {
            if let Err(error) = self.source.block(output) {
                self.poisoned = true;
                return Err(error);
            }
        } else {
            let source_capacity = self.source.max_block_len();
            for value in output {
                // A source admitted with --block 1 is still legal. These splits
                // change neither source event times nor mechanical arithmetic.
                for chunk in self.input.chunks_mut(source_capacity) {
                    if let Err(error) = self.source.block(chunk) {
                        self.poisoned = true;
                        return Err(error);
                    }
                }
                match self.filter.preview(&self.input) {
                    Ok(frame) => *value = frame[0],
                    Err(what) => {
                        self.poisoned = true;
                        return Err(RenderError::Control { what });
                    }
                }
                self.filter.commit();
            }
        }
        self.completed += samples;
        Ok(())
    }
}
impl<S: PressureRenderer> PressureRenderer for DecimatedRenderer<S> {
    fn samples_rendered(&self) -> u64 { Self::samples_rendered(self) }
    fn max_block_len(&self) -> usize { Self::max_block_len(self) }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> { Self::validate_sample_rate(self, rate) }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> { Self::validate_sample_count(self, samples) }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> { Self::block(self, output) }
}
fn sizing(what: &'static str) -> RenderError { RenderError::Sizing { what } }
