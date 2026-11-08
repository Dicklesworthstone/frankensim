//! Retained BEM microphone/beamformer waveforms through the existing WAV owner.

use crate::pcm_wav::{MAX_WAV_CHANNELS, observation::PressureFrameRenderer};
use crate::render::RenderError;
use fs_bem::observation_transfers::{
    HarmonicPressureChannel, HarmonicPressureSampleLimits, HarmonicSampleClock,
    PreparedHarmonicPressureSignal, TransferError,
};

/// Stream one original BEM realization's directional/delayed pressure channels.
///
/// Fixed source draws, native angular frequencies, original p/u microphone
/// projection and signed causal taps belong to the borrowed signal. This
/// adapter adds one bounded reusable callback buffer, not another propagator,
/// source envelope, gain normalization or PCM encoder.
///
/// Callbacks publish all channels together. A refusal or cancellation leaves
/// both caller output and frame counter unchanged, permitting a safe retry.
/// Successful callbacks allocate nothing and preserve block partition bits.
pub struct BemHarmonicPressureFrameRenderer<'a> {
    signal: &'a PreparedHarmonicPressureSignal<'a>,
    clock: HarmonicSampleClock,
    limits: HarmonicPressureSampleLimits,
    scratch: Vec<f64>,
    completed: u64,
}

impl<'a> BemHarmonicPressureFrameRenderer<'a> {
    pub fn new(
        signal: &'a PreparedHarmonicPressureSignal<'a>,
        clock: HarmonicSampleClock,
        limits: HarmonicPressureSampleLimits,
    ) -> Result<Self, RenderError> {
        let channels = signal.channel_count();
        if channels == 0 || channels > usize::from(MAX_WAV_CHANNELS) || limits.maximum_frames == 0 {
            return Err(sizing(
                "BEM harmonic output needs bounded channels and callback frames",
            ));
        }
        signal
            .validate_sample_block(clock, 0, limits.maximum_frames, limits)
            .map_err(RenderError::BemObservation)?;
        let count = channels
            .checked_mul(limits.maximum_frames)
            .ok_or_else(|| sizing("BEM harmonic callback storage count overflow"))?;
        let mut scratch = Vec::new();
        scratch
            .try_reserve_exact(count)
            .map_err(|_| sizing("BEM harmonic callback allocation failed"))?;
        scratch.resize(count, 0.0);
        Ok(Self {
            signal,
            clock,
            limits,
            scratch,
            completed: 0,
        })
    }

    pub fn signal(&self) -> &PreparedHarmonicPressureSignal<'a> {
        self.signal
    }

    /// Preserve this authored channel order beside exported WAV bytes.
    pub fn channels(&self) -> &[HarmonicPressureChannel] {
        self.signal.pressure_array().channels()
    }

    pub const fn clock(&self) -> HarmonicSampleClock {
        self.clock
    }

    pub fn frame_block_cancellable(
        &mut self,
        output: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), RenderError> {
        poll(&mut cancelled)?;
        if output.is_empty() {
            return Err(RenderError::EmptyBlock);
        }
        let channels = self.signal.channel_count();
        if !output.len().is_multiple_of(channels) || output.len() > self.scratch.len() {
            return Err(sizing(
                "BEM harmonic output needs complete frames within callback capacity",
            ));
        }
        let frames = u64::try_from(output.len() / channels)
            .map_err(|_| sizing("BEM harmonic frame count exceeds u64"))?;
        self.signal
            .sample_pressure_into_cancellable(
                self.clock,
                self.completed,
                &mut self.scratch[..output.len()],
                self.limits,
                &mut cancelled,
            )
            .map_err(RenderError::BemObservation)?;
        poll(&mut cancelled)?;
        output.copy_from_slice(&self.scratch[..output.len()]);
        self.completed += frames;
        Ok(())
    }
}

impl PressureFrameRenderer for BemHarmonicPressureFrameRenderer<'_> {
    fn channel_count(&self) -> usize {
        self.signal.channel_count()
    }

    fn frames_rendered(&self) -> u64 {
        self.completed
    }

    fn max_frame_block_len(&self) -> usize {
        self.limits.maximum_frames
    }

    fn validate_frame_rate(&self, rate: u32) -> Result<(), RenderError> {
        if rate != self.clock.sample_rate_hz() {
            return Err(sizing("BEM harmonic source and WAV sample clocks differ"));
        }
        Ok(())
    }

    fn validate_frame_count(&self, frames: u64) -> Result<(), RenderError> {
        self.signal
            .realization()
            .validate_sample_range(self.clock, self.completed, frames)
            .map_err(RenderError::BemObservation)
    }

    fn frame_block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        self.frame_block_cancellable(output, || false)
    }
}

fn sizing(what: &'static str) -> RenderError {
    RenderError::Sizing { what }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RenderError> {
    if cancelled() {
        Err(RenderError::BemObservation(TransferError::Cancelled))
    } else {
        Ok(())
    }
}
