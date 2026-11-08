//! Sampled directional/delayed pressure channels of one fixed physical draw.
//!
//! Projection and delays are prepared by the existing pressure-array owner.
//! No new factor draw, pressure-to-velocity approximation or source onset is
//! introduced when a harmonic field becomes an audible pressure signal.

use super::super::realization::phase;
use super::*;

/// Bounded sampling of processed channels, not new spatial receivers.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicPressureSampleLimits {
    pub maximum_frames: usize,
    /// Frames times output channels times retained physical frequencies.
    pub maximum_channel_line_samples: usize,
}

/// Pressure-equivalent microphone/beamformer waveforms of one realization.
///
/// Every channel uses the SAME fixed real/proper factor draws already present
/// in the borrowed realization. The original array owner prepares p/u
/// projection, signed receiver/tap summation and exp(+i omega delay) exactly
/// once in DD. Thus nulls, reactive velocity and fractional delays survive
/// before the final scalar sample, without resampling or phase interpolation.
///
/// The retained pressure array has only a deterministic mean bank; its window
/// moments describe this fixed waveform, not the source ensemble's moments.
/// Original receivers and source identity remain distinct from output channels.
pub struct PreparedHarmonicPressureSignal<'a> {
    realization: &'a PreparedHarmonicRealization<'a>,
    array: PreparedHarmonicPressureArray<'a>,
}

impl<'a> PreparedHarmonicPressureSignal<'a> {
    /// Prepare authored taps. None observes physical pressure; Some selects one
    /// existing ideal FirstOrder pattern per ORIGINAL physical receiver.
    pub fn new(
        realization: &'a PreparedHarmonicRealization<'a>,
        channels: &[HarmonicPressureChannel],
        patterns: Option<&[FirstOrder]>,
        limits: HarmonicPressureArrayLimits,
    ) -> Result<Self, TransferError> {
        Self::new_cancellable(realization, channels, patterns, limits, || false)
    }

    /// Uses the same channel, pattern, work and storage admission as windowed
    /// pressure arrays. All preparation is atomic and leaves its source intact.
    pub fn new_cancellable(
        realization: &'a PreparedHarmonicRealization<'a>,
        channels: &[HarmonicPressureChannel],
        patterns: Option<&[FirstOrder]>,
        limits: HarmonicPressureArrayLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        let array = PreparedHarmonicPressureArray::prepare(
            realization.pressure_input(),
            channels,
            patterns,
            limits,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(Self { realization, array })
    }

    /// Original physical waveform and its unchanged shared factor draws.
    pub fn realization(&self) -> &PreparedHarmonicRealization<'a> {
        self.realization
    }

    /// Original channel definitions, patterns, scope and exact window moments
    /// of this deterministic waveform. No ensemble variance is inferred here.
    pub fn pressure_array(&self) -> &PreparedHarmonicPressureArray<'a> {
        &self.array
    }

    pub fn channel_count(&self) -> usize {
        self.array.channels.len()
    }

    pub fn line_count(&self) -> usize {
        self.realization.line_count()
    }

    /// Preflight a callback before allocating a buffer or advancing a writer.
    /// The complete stream horizon can separately be checked through the
    /// original realization's validate_sample_range on the same clock.
    pub fn validate_sample_block(
        &self,
        clock: HarmonicSampleClock,
        first_sample: u64,
        frames: usize,
        limits: HarmonicPressureSampleLimits,
    ) -> Result<(), TransferError> {
        self.admit_block(clock, first_sample, frames, limits, &mut || false)
    }

    fn admit_block(
        &self,
        clock: HarmonicSampleClock,
        first_sample: u64,
        frames: usize,
        limits: HarmonicPressureSampleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        poll(cancelled)?;
        admit(frames, limits.maximum_frames)?;
        let work = frames
            .checked_mul(self.channel_count())
            .and_then(|n| n.checked_mul(self.line_count()))
            .ok_or_else(overflow)?;
        admit(work, MAX_PRODUCTS.min(limits.maximum_channel_line_samples))?;
        let count = u64::try_from(frames).map_err(|_| overflow())?;
        self.realization
            .validate_sample_range_cancellable(clock, first_sample, count, cancelled)
    }

    /// Frame-interleaved physical Pa-equivalent output in authored channel
    /// order. Every original frequency must be below Nyquist, including zero
    /// taps and coherent nulls. Empty channels admit only an empty buffer.
    ///
    /// Sampling allocates nothing. Absolute integer indices share the source
    /// clock and are independent of block partition. Admission failures leave
    /// output untouched; discard it after runtime overflow or cancellation,
    /// which may have written a prefix.
    pub fn sample_pressure_into(
        &self,
        clock: HarmonicSampleClock,
        first_sample: u64,
        output: &mut [f64],
        limits: HarmonicPressureSampleLimits,
    ) -> Result<(), TransferError> {
        self.sample_pressure_into_cancellable(clock, first_sample, output, limits, || false)
    }

    pub fn sample_pressure_into_cancellable(
        &self,
        clock: HarmonicSampleClock,
        first_sample: u64,
        output: &mut [f64],
        limits: HarmonicPressureSampleLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        poll(&mut cancelled)?;
        let channels = self.channel_count();
        if channels == 0 {
            if !output.is_empty() {
                return Err(bad("empty harmonic pressure channels require empty output"));
            }
            return self.admit_block(clock, first_sample, 0, limits, &mut cancelled);
        }
        if !output.len().is_multiple_of(channels) {
            return Err(bad(
                "harmonic pressure output needs complete channel frames",
            ));
        }
        self.admit_block(
            clock,
            first_sample,
            output.len() / channels,
            limits,
            &mut cancelled,
        )?;
        let mut sums = [Dd::ZERO; MAX_PRESSURE_ARRAY_CHANNELS];
        for (offset, frame) in output.chunks_exact_mut(channels).enumerate() {
            poll(&mut cancelled)?;
            let time = clock.time_at(first_sample + offset as u64)?;
            sums[..channels].fill(Dd::ZERO);
            // This owner is prepared from one deterministic field bank.
            for line in &self.array.banks[0].lines {
                poll(&mut cancelled)?;
                let [real, imaginary] = phase(line.omega, time)?;
                for (sum, amplitude) in sums[..channels].iter_mut().zip(&line.amplitudes) {
                    *sum = *sum + amplitude[0] * real - amplitude[1] * imaginary;
                }
            }
            for (sample, value) in frame.iter_mut().zip(&sums[..channels]) {
                if !finite(*value) || !value.to_f64().is_finite() {
                    return Err(bad("harmonic pressure signal sample is unrepresentable"));
                }
                *sample = value.to_f64();
            }
        }
        poll(&mut cancelled)
    }
}
