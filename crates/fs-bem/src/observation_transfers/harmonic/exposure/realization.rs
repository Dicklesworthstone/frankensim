//! Physical waveforms from the original retained BEM pressure/velocity fields.
//!
//! A supplied draw multiplies its whole factor across all frequencies and
//! receivers. Draws never change during sampling. No random distribution,
//! propagation delay, transient onset or plane-wave velocity is invented.

use super::*;
use crate::near_field::MAX_RECEIVERS;
use fs_math::det;

/// One fixed caller-supplied coordinate of the authored harmonic ensemble.
///
/// Across an ensemble, real coordinates have E[xi]=0 and E[xi²]=1; proper
/// coordinates have E[z]=E[z²]=0 and E[|z|²]=1. A single draw need not have
/// unit magnitude. Values multiply the peak phasors directly, without a
/// conjugation, normalization, sqrt(2) or extra peak-to-RMS conversion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HarmonicFactorDraw {
    Real(f64),
    ProperComplex(C64),
}

impl HarmonicFactorDraw {
    fn parts(self, law: HarmonicExposureFactorStatistics) -> Result<[Dd; 2], TransferError> {
        let value = match (self, law) {
            (Self::Real(value), HarmonicExposureFactorStatistics::Real) => C64::from_re(value),
            (Self::ProperComplex(value), HarmonicExposureFactorStatistics::ProperComplex) => value,
            _ => return Err(bad("BEM realization draw does not match its factor law")),
        };
        if !value.re.is_finite() || !value.im.is_finite() {
            return Err(bad("BEM realization requires finite factor draws"));
        }
        Ok([Dd::from_f64(value.re), Dd::from_f64(value.im)])
    }
}

/// Bounds on preparing one coherent member of an existing exposure ensemble.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicRealizationLimits {
    /// Sum of retained mean/factor lines BEFORE grouping across factors.
    pub maximum_lines: usize,
    /// Four complex field products per receiver and supplied line.
    pub maximum_field_products: usize,
    /// Peak logical heap payload: fields, line headers, sort entries and draws.
    /// Borrowed owners, allocator overhead and fixed-size stack are excluded.
    pub maximum_workspace_bytes: usize,
}

/// A uniform clock on the same physical epoch as the source phasors.
/// Samples occur at `epoch_s + index / sample_rate_hz`; propagation is already
/// in the retained fields. DD preserves short intervals at large epochs.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicSampleClock {
    epoch_s: Dd,
    sample_rate_hz: u32,
}

impl HarmonicSampleClock {
    pub fn new(epoch_s: Dd, sample_rate_hz: u32) -> Result<Self, TransferError> {
        if sample_rate_hz == 0 || !finite(epoch_s) {
            return Err(bad(
                "BEM sample clock needs a finite epoch and positive rate",
            ));
        }
        let clock = Self {
            epoch_s: Dd::from_pair(epoch_s.hi, epoch_s.lo),
            sample_rate_hz,
        };
        clock.validate_index(0)?;
        Ok(clock)
    }

    pub const fn epoch_s(self) -> Dd {
        self.epoch_s
    }

    pub const fn sample_rate_hz(self) -> u32 {
        self.sample_rate_hz
    }

    /// Retain the complete integer counter, including indices beyond 2^53.
    pub fn time_at(self, index: u64) -> Result<Dd, TransferError> {
        let integer = Dd::from_f64((index >> 32) as f64) * Dd::from_f64(4_294_967_296.0)
            + Dd::from_f64((index & 0xffff_ffff) as f64);
        let time = self.epoch_s + integer / Dd::from_f64(f64::from(self.sample_rate_hz));
        if !finite(time) {
            return Err(bad("BEM sample time is unrepresentable"));
        }
        Ok(time)
    }

    fn validate_index(self, index: u64) -> Result<(), TransferError> {
        let time = self.time_at(index)?;
        let next = time + Dd::ONE / Dd::from_f64(f64::from(self.sample_rate_hz));
        if !finite(next) || !time.lt(next) {
            return Err(bad("BEM sample interval is unresolved on the source clock"));
        }
        Ok(())
    }
}

/// Work for one caller-owned, sample-major output buffer.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicSampleLimits {
    pub maximum_frames: usize,
    /// Frames times receivers times distinct realized physical frequencies.
    /// Each term updates four real pressure/velocity components.
    pub maximum_receiver_line_samples: usize,
}

/// Simultaneous instantaneous pressure and actual Cartesian particle velocity.
/// Values are physical point samples, not RMS amplitudes or window averages.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HarmonicFieldSample {
    pub pressure_pa: f64,
    pub particle_velocity_m_s: [f64; 3],
}

/// One coherent waveform realization of the original BEM source ensemble.
///
/// Preparation contracts fixed factor draws once, retaining DD coefficients
/// and exactly grouping equal native angular frequencies. Close frequencies
/// remain distinct. Even zero draws and exact coherent nulls retain their
/// physical lines for subsequent Nyquist admission. Sampling allocates nothing
/// and evaluates each temporal phase once per frequency across all receivers.
/// The exposure may be dropped; original physical identity remains borrowed.
pub struct PreparedHarmonicRealization<'a> {
    reference: &'a PreparedHarmonicObservation,
    source_clock: &'a str,
    draws: Vec<HarmonicFactorDraw>,
    lines: Vec<Line>,
    workspace_bytes: usize,
}

struct Input<'a> {
    line: &'a Line,
    draw: [Dd; 2],
    order: usize,
}

impl<'a> PreparedHarmonicExposure<'a> {
    /// Apply one fixed draw per original factor, in its declared order.
    pub fn realize(
        &self,
        draws: &[HarmonicFactorDraw],
        limits: HarmonicRealizationLimits,
    ) -> Result<PreparedHarmonicRealization<'a>, TransferError> {
        self.realize_cancellable(draws, limits, || false)
    }

    /// Bounded preparation; refusal/cancellation returns no partial owner.
    pub fn realize_cancellable(
        &self,
        draws: &[HarmonicFactorDraw],
        limits: HarmonicRealizationLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedHarmonicRealization<'a>, TransferError> {
        poll(&mut cancelled)?;
        if draws.len() != self.factors.len() {
            return Err(bad("BEM realization needs exactly one draw per factor"));
        }
        let mut count = self.mean.lines.len();
        for (&draw, bank) in draws.iter().zip(&self.factors) {
            poll(&mut cancelled)?;
            draw.parts(bank.statistics)?;
            count = count.checked_add(bank.lines.len()).ok_or_else(overflow)?;
        }
        let receivers = self.reference.receivers.len();
        admit(count, MAX_EXPOSURE_LINES.min(limits.maximum_lines))?;
        let products = count
            .checked_mul(receivers)
            .and_then(|n| n.checked_mul(4))
            .ok_or_else(overflow)?;
        admit(products, MAX_PRODUCTS.min(limits.maximum_field_products))?;
        let per_line = receivers
            .checked_mul(core::mem::size_of::<[Dd; 8]>())
            .and_then(|n| n.checked_add(core::mem::size_of::<Line>()))
            .and_then(|n| n.checked_add(core::mem::size_of::<Input<'_>>()))
            .ok_or_else(overflow)?;
        let workspace_bytes = count
            .checked_mul(per_line)
            .and_then(|n| {
                n.checked_add(
                    draws
                        .len()
                        .checked_mul(core::mem::size_of::<HarmonicFactorDraw>())?,
                )
            })
            .ok_or_else(overflow)?;
        admit(workspace_bytes, limits.maximum_workspace_bytes)?;

        let mut inputs = storage(count)?;
        for (index, bank) in core::iter::once(&self.mean)
            .chain(&self.factors)
            .enumerate()
        {
            let draw = if index == 0 {
                [Dd::ONE, Dd::ZERO]
            } else {
                draws[index - 1].parts(bank.statistics)?
            };
            for line in &bank.lines {
                poll(&mut cancelled)?;
                inputs.push(Input {
                    line,
                    draw,
                    order: inputs.len(),
                });
            }
        }
        // Equal-frequency addition preserves mean, then original factor order.
        inputs.sort_unstable_by(|a, b| {
            a.line
                .omega
                .total_cmp(&b.line.omega)
                .then_with(|| a.order.cmp(&b.order))
        });
        let mut lines: Vec<Line> = storage(count)?;
        for input in inputs {
            poll(&mut cancelled)?;
            if lines
                .last()
                .is_none_or(|last| last.omega.to_bits() != input.line.omega.to_bits())
            {
                let mut fields = storage(receivers)?;
                fields.resize(receivers, [Dd::ZERO; 8]);
                lines.push(Line {
                    omega: input.line.omega,
                    fields,
                });
            }
            let fields = &mut lines.last_mut().expect("frequency group exists").fields;
            let [real, imaginary] = input.draw;
            for (out, field) in fields.iter_mut().zip(&input.line.fields) {
                for c in 0..4 {
                    out[2 * c] = out[2 * c] + field[2 * c] * real - field[2 * c + 1] * imaginary;
                    out[2 * c + 1] =
                        out[2 * c + 1] + field[2 * c] * imaginary + field[2 * c + 1] * real;
                }
            }
        }
        if lines
            .iter()
            .flat_map(|line| &line.fields)
            .flatten()
            .any(|v| !finite(*v))
        {
            return Err(bad("BEM realized field is unrepresentable"));
        }
        let mut retained_draws = storage(draws.len())?;
        retained_draws.extend_from_slice(draws);
        poll(&mut cancelled)?;
        Ok(PreparedHarmonicRealization {
            reference: self.reference,
            source_clock: self.source_clock,
            draws: retained_draws,
            lines,
            workspace_bytes,
        })
    }
}

impl PreparedHarmonicRealization<'_> {
    pub(super) fn pressure_input(&self) -> super::array::PressureInput<'_, '_> {
        super::array::PressureInput {
            reference: self.reference,
            source_clock: self.source_clock,
            mean: &self.lines,
            factors: &[],
        }
    }

    pub fn receivers(&self) -> &[ReceiverSpec] {
        &self.reference.receivers
    }

    pub fn sources(&self) -> &[SourceId] {
        &self.reference.sources
    }

    pub fn source_clock(&self) -> &str {
        self.source_clock
    }

    pub fn producer_hash(&self) -> ContentHash {
        self.reference.producer_hash
    }

    pub fn producer_uncertainty(&self) -> ProducerUncertainty {
        ProducerUncertainty::Unresolved {
            source_hash: self.producer_hash(),
        }
    }

    pub fn factor_draws(&self) -> &[HarmonicFactorDraw] {
        &self.draws
    }

    pub fn line_count(&self) -> usize {
        self.lines.len()
    }

    /// Admitted peak logical heap payload of preparation.
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Admit a whole range before a streaming consumer changes its output.
    /// Every original physical line must be strictly below Nyquist, including
    /// zero-draw and coherently cancelled lines. Empty ranges are checked too.
    pub fn validate_sample_range(
        &self,
        clock: HarmonicSampleClock,
        first_sample: u64,
        sample_count: u64,
    ) -> Result<(), TransferError> {
        self.validate_sample_range_cancellable(clock, first_sample, sample_count, || false)
    }

    pub fn validate_sample_range_cancellable(
        &self,
        clock: HarmonicSampleClock,
        first_sample: u64,
        sample_count: u64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        poll(&mut cancelled)?;
        let end = first_sample
            .checked_add(sample_count)
            .ok_or_else(|| bad("BEM sample counter overflow"))?;
        clock.validate_index(first_sample)?;
        clock.validate_index(end)?;
        let first_time = clock.time_at(first_sample)?;
        let end_time = clock.time_at(end)?;
        let nyquist =
            fs_fft::harmonic_window::angular_frequency(0.5 * f64::from(clock.sample_rate_hz))
                .map_err(window_error)?;
        for line in &self.lines {
            poll(&mut cancelled)?;
            if !Dd::from_f64(line.omega).lt(nyquist) {
                return Err(bad("BEM harmonic line is at or above output Nyquist"));
            }
            phase(line.omega, first_time)?;
            phase(line.omega, end_time)?;
        }
        poll(&mut cancelled)
    }

    /// Sample all receivers in `output[frame * receiver_count + receiver]`.
    /// Shape, range, work and bandwidth are checked before any write. Runtime
    /// overflow or cancellation may leave a prefix; discard output on error.
    /// Absolute integer indices make arbitrary block partitioning bit-identical
    /// and avoid oscillator drift or new draws at callback boundaries.
    pub fn sample_fields_into(
        &self,
        clock: HarmonicSampleClock,
        first_sample: u64,
        output: &mut [HarmonicFieldSample],
        limits: HarmonicSampleLimits,
    ) -> Result<(), TransferError> {
        self.sample_fields_into_cancellable(clock, first_sample, output, limits, || false)
    }

    /// Allocation-free p/u synthesis; one DD phase per line serves the array.
    pub fn sample_fields_into_cancellable(
        &self,
        clock: HarmonicSampleClock,
        first_sample: u64,
        output: &mut [HarmonicFieldSample],
        limits: HarmonicSampleLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), TransferError> {
        poll(&mut cancelled)?;
        let receivers = self.reference.receivers.len();
        if !output.len().is_multiple_of(receivers) {
            return Err(bad("BEM sample output must contain whole receiver frames"));
        }
        let frames = output.len() / receivers;
        admit(frames, limits.maximum_frames)?;
        let work = output
            .len()
            .checked_mul(self.lines.len())
            .ok_or_else(overflow)?;
        admit(work, MAX_PRODUCTS.min(limits.maximum_receiver_line_samples))?;
        let count = u64::try_from(frames).map_err(|_| overflow())?;
        self.validate_sample_range_cancellable(clock, first_sample, count, &mut cancelled)?;

        let mut sums = [[Dd::ZERO; 4]; MAX_RECEIVERS];
        for (offset, frame) in output.chunks_exact_mut(receivers).enumerate() {
            poll(&mut cancelled)?;
            let time = clock.time_at(first_sample + offset as u64)?;
            sums[..receivers].fill([Dd::ZERO; 4]);
            for line in &self.lines {
                poll(&mut cancelled)?;
                let [real, imaginary] = phase(line.omega, time)?;
                for (sum, field) in sums[..receivers].iter_mut().zip(&line.fields) {
                    for c in 0..4 {
                        sum[c] = sum[c] + field[2 * c] * real - field[2 * c + 1] * imaginary;
                    }
                }
            }
            for (sample, sum) in frame.iter_mut().zip(&sums[..receivers]) {
                let value = sum.map(Dd::to_f64);
                if sum.iter().any(|v| !finite(*v)) || value.iter().any(|v| !v.is_finite()) {
                    return Err(bad("BEM sampled field is unrepresentable"));
                }
                *sample = HarmonicFieldSample {
                    pressure_pa: value[0],
                    particle_velocity_m_s: [value[1], value[2], value[3]],
                };
            }
        }
        poll(&mut cancelled)
    }
}

pub(super) fn phase(omega: f64, time: Dd) -> Result<[Dd; 2], TransferError> {
    // Keep the original rad/s clock, never a rounded round-trip through Hz.
    // This is the same DD transcendental owner as fs-fft harmonic windows.
    let angle = Dd::from_f64(omega) * time;
    if !finite(angle) {
        return Err(bad("BEM sampled harmonic phase is unrepresentable"));
    }
    let (sine, cosine) = det::sin_cos_dd(angle);
    if !finite(sine) || !finite(cosine) {
        return Err(bad("BEM sampled harmonic phase is unrepresentable"));
    }
    Ok([cosine, -sine])
}
