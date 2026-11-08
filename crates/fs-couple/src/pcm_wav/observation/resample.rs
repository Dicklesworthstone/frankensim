//! Causal, bounded fractional-rate observation of simultaneous interval means.
//!
//! An input frame is the ORIGINAL mean on `[i*dt, (i+1)*dt]`, located at its
//! center and available at its end. Output centers are `(j+1/2)/output_rate`.
//! Only observation is filtered: the source still advances through every
//! interval covering each original output window, including the final one.
//! A symmetric finite sinc with half-support `H` observes the input sequence
//! with delay `(H+1/2)*dt`; negative source indices are zero. No future source
//! frames, reflected edges, padded physical tail or delay removal are used.
//!
//! The four-term Blackman-Harris window uses the same family as `Decimator`,
//! with its tiny endpoint pedestal removed so phase wrap is continuous.
//! `H = ceil(4 / (dt * (stopband-passband)))`; the cutoff is the transition's
//! midpoint. Each phase row has unity DC gain. Adjacent cold-prepared rows are
//! interpolated linearly, sharing one phase across all channels. Preparation
//! owns all allocation and transcendental evaluation; push/render do neither.
//!
//! This is a finite-FIR approximation, not a certified brick-wall filter.
//! Transition bandwidth, phase-table spacing and retained latency are explicit.
//! It reconstructs the sequence of interval means: it does NOT divide out their
//! physical averaging sinc, restore content already aliased by the original
//! time step, certify source convergence or filter forces/material laws.
//!
//! Reference: J. O. Smith, *Physical Audio Signal Processing*,
//! [Windowed Sinc Interpolation](https://www.dsprelated.com/freebooks/pasp/Windowed_Sinc_Interpolation.html).

use fs_math::det;
use std::f64::consts::PI;

const EXACT_INTEGER_LIMIT: u64 = 1 << 53;
// Limits complement caller budgets, including deliberately enormous inputs.
const MAX_RETAINED_VALUES: usize = 1 << 23;
const MAX_PREPARATION_TERMS: usize = 1 << 29;
const MAX_OUTPUT_TERMS: usize = 1 << 22;

/// Explicit bandwidth, fractional-phase resolution and cold/warm work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BandlimitedObservationPolicy {
    /// End of the retained passband [Hz], strictly positive.
    pub passband_hz: f64,
    /// Start of the stopband [Hz], at most either clock's Nyquist frequency.
    pub stopband_hz: f64,
    /// Number of fractional intervals in [0,1]; 16..=16384. The table retains
    /// one extra endpoint row. Larger values reduce interpolation error, not
    /// the finite window's stopband leakage.
    pub phase_count: usize,
    /// Budget `taps * channels * (phase_count+1)`. Coefficients are physically
    /// shared across channels; admission charges their complete channel work.
    pub max_preparation_terms: usize,
    /// Retained scalar input-history values plus one transactional output frame.
    pub max_history_values: usize,
    /// Budget `taps * channels` per emitted frame, independent of preparation.
    pub max_output_terms: usize,
}

/// Actual admitted filter dimensions and timing, not an attenuation certificate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BandlimitedObservationReport {
    pub input_period_s: f64,
    pub output_rate_hz: u32,
    pub channels: usize,
    pub passband_hz: f64,
    pub stopband_hz: f64,
    pub taps: usize,
    pub half_support: usize,
    pub phase_count: usize,
    /// Spacing of the table in ORIGINAL input-frame units.
    pub phase_spacing_input_samples: f64,
    pub latency_s: f64,
    pub latency_output_samples: f64,
    /// Includes extra accepted source frames between an output's center/end.
    pub history_frames: usize,
    /// History plus the one-frame output scratch allocation, in scalar values.
    pub history_values: usize,
    pub preparation_terms: usize,
    pub output_terms: usize,
}

/// One finite observation window, starting with zero prehistory.
///
/// For each output, push through `required_source_intervals()`, then call
/// `render_next`. The required count is absolute and can be unchanged between
/// outputs during upsampling. It always covers the original physical output
/// interval, even when startup yields zero. Feeding ahead of that boundary is
/// refused so retained memory is independent of the total duration.
#[derive(Clone, Debug)]
pub struct BandlimitedResampler {
    policy: BandlimitedObservationPolicy,
    report: BandlimitedObservationReport,
    physical_width: f64,
    coverage_width: f64,
    output_samples: u64,
    source_intervals: usize,
    accepted: usize,
    rendered: u64,
    coefficients: Vec<f64>,
    history: Vec<f64>,
    output: Vec<f64>,
}

impl BandlimitedResampler {
    /// Prepare all coefficients/history without advancing the physical source.
    ///
    /// Refuses invalid clocks/bands, inadequate budgets and a horizon whose
    /// floating clock cannot resolve the requested fractional-phase spacing.
    /// The latter is normally far beyond the caller's physical-work budget.
    pub fn new(
        input_period_s: f64,
        output_rate_hz: u32,
        channels: usize,
        output_samples: u64,
        policy: BandlimitedObservationPolicy,
    ) -> Result<Self, &'static str> {
        Self::new_cancellable(
            input_period_s,
            output_rate_hz,
            channels,
            output_samples,
            policy,
            &mut || false,
        )
    }

    /// Same cold preparation with cancellation checked before allocation and
    /// at most every 256 coefficient evaluations. No physical source is owned
    /// or advanced by either construction path.
    pub fn new_cancellable(
        input_period_s: f64,
        output_rate_hz: u32,
        channels: usize,
        output_samples: u64,
        policy: BandlimitedObservationPolicy,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, &'static str> {
        if cancelled() {
            return Err("bandlimited observation preparation cancelled");
        }
        let physical_width = input_period_s * f64::from(output_rate_hz);
        if !input_period_s.is_finite()
            || input_period_s <= 0.0
            || output_rate_hz == 0
            || !physical_width.is_finite()
            || physical_width <= 0.0
            || channels == 0
            || output_samples == 0
            || output_samples > EXACT_INTEGER_LIMIT / 2
        {
            return Err(
                "bandlimited observation needs positive finite clocks, channels and representable output centers",
            );
        }
        let input_rate = input_period_s.recip();
        let nyquist = 0.5 * input_rate.min(f64::from(output_rate_hz));
        if !input_rate.is_finite()
            || !policy.passband_hz.is_finite()
            || !policy.stopband_hz.is_finite()
            || policy.passband_hz <= 0.0
            || policy.passband_hz >= policy.stopband_hz
            || policy.stopband_hz > nyquist
            || !(16..=16384).contains(&policy.phase_count)
        {
            return Err(
                "bandlimited observation requires 0 < passband < stopband <= both Nyquists and 16..=16384 phases",
            );
        }
        let half = (4.0 / (input_period_s * (policy.stopband_hz - policy.passband_hz))).ceil();
        if !half.is_finite() || half < 1.0 || half > MAX_RETAINED_VALUES as f64 {
            return Err("bandlimited observation transition requires excessive finite support");
        }
        let half_support = half as usize;
        let taps = 2 * half_support + 1;
        let rows = policy.phase_count + 1;
        let output_terms = taps
            .checked_mul(channels)
            .ok_or("bandlimited observation output work overflows")?;
        let preparation_terms = output_terms
            .checked_mul(rows)
            .ok_or("bandlimited observation preparation work overflows")?;
        let table_values = taps
            .checked_mul(rows)
            .ok_or("bandlimited observation phase table overflows")?;
        // All original source intervals covering the output END are accepted
        // before observation at its CENTER. Keep those in addition to the
        // delayed 2H+1 stencil. The small margin also covers clock roundoff.
        let ahead = (0.5 / physical_width).ceil();
        if !ahead.is_finite() || ahead > MAX_RETAINED_VALUES as f64 {
            return Err("bandlimited observation source/output ratio exceeds retained history");
        }
        let history_frames = taps + ahead as usize + 3;
        let history_values = (history_frames + 1)
            .checked_mul(channels)
            .ok_or("bandlimited observation history size overflows")?;
        if output_terms > policy.max_output_terms
            || output_terms > MAX_OUTPUT_TERMS
            || preparation_terms > policy.max_preparation_terms
            || preparation_terms > MAX_PREPARATION_TERMS
            || table_values > MAX_RETAINED_VALUES
            || history_values > policy.max_history_values
            || history_values > MAX_RETAINED_VALUES
        {
            return Err(
                "bandlimited observation exceeds preparation, history or per-output work budget",
            );
        }
        // Match the existing NSF same-clock admission ONLY for physical-window
        // coverage. The actual physical dt is retained everywhere in the FIR.
        let output_period_s = 1.0 / f64::from(output_rate_hz);
        let matched = (input_period_s - output_period_s).abs()
            <= 8.0 * f64::EPSILON * input_period_s.abs().max(output_period_s.abs());
        let coverage_width = if matched { 1.0 } else { physical_width };
        let source_intervals = coverage(output_samples, coverage_width)?;
        let end_coordinate = output_samples as f64 / physical_width;
        // Eight roundoffs conservatively cover index-center division, width
        // multiplication and support subtraction. Keep clock granularity well
        // below table resolution, rather than pretending huge indices are exact.
        if !end_coordinate.is_finite()
            || 8.0 * f64::EPSILON * (end_coordinate + half + 1.0)
                > 1.0 / (16.0 * policy.phase_count as f64)
        {
            return Err(
                "bandlimited observation horizon cannot resolve its fractional-phase table",
            );
        }
        let latency_s = (half + 0.5) * input_period_s;
        let latency_output_samples = (half + 0.5) * physical_width;
        if !latency_s.is_finite() || !latency_output_samples.is_finite() {
            return Err("bandlimited observation latency is not finite");
        }
        let report = BandlimitedObservationReport {
            input_period_s,
            output_rate_hz,
            channels,
            passband_hz: policy.passband_hz,
            stopband_hz: policy.stopband_hz,
            taps,
            half_support,
            phase_count: policy.phase_count,
            phase_spacing_input_samples: 1.0 / policy.phase_count as f64,
            latency_s,
            latency_output_samples,
            history_frames,
            history_values,
            preparation_terms,
            output_terms,
        };
        if cancelled() {
            return Err("bandlimited observation preparation cancelled");
        }
        let mut coefficients = zeros(table_values)?;
        let cutoff = 0.5 * (policy.passband_hz + policy.stopband_hz) * input_period_s;
        for phase in 0..policy.phase_count {
            let fraction = phase as f64 / policy.phase_count as f64;
            let row = &mut coefficients[phase * taps..(phase + 1) * taps];
            let mut sum = 0.0;
            for (tap, value) in row.iter_mut().enumerate() {
                if tap % 256 == 0 && cancelled() {
                    return Err("bandlimited observation preparation cancelled");
                }
                *value = kernel(fraction - (tap as f64 - half), half, cutoff);
                sum += *value;
            }
            if !sum.is_finite() || sum <= 0.0 {
                return Err("bandlimited observation phase has nonpositive normalization");
            }
            for value in row {
                *value /= sum;
            }
        }
        // At phase=1, the identical physical stencil is phase=0 with its
        // anchor shifted one source frame. Copy it exactly (both endpoints of
        // the window are zero), instead of independently rounding the seam.
        for tap in 1..taps {
            coefficients[policy.phase_count * taps + tap] = coefficients[tap - 1];
        }
        if cancelled() {
            return Err("bandlimited observation preparation cancelled");
        }
        Ok(Self {
            policy,
            report,
            physical_width,
            coverage_width,
            output_samples,
            source_intervals,
            accepted: 0,
            rendered: 0,
            coefficients,
            history: zeros(history_frames * channels)?,
            output: zeros(channels)?,
        })
    }

    pub const fn policy(&self) -> BandlimitedObservationPolicy {
        self.policy
    }
    pub const fn report(&self) -> BandlimitedObservationReport {
        self.report
    }
    pub const fn accepted_source_intervals(&self) -> usize {
        self.accepted
    }
    pub const fn rendered_samples(&self) -> u64 {
        self.rendered
    }
    /// Total physical frames, including final partial output-window coverage.
    pub const fn source_intervals(&self) -> usize {
        self.source_intervals
    }
    pub const fn output_samples(&self) -> u64 {
        self.output_samples
    }

    /// Absolute physical interval count needed before the NEXT output frame.
    /// Exhaustion is a refusal; it never requests a padded tail.
    pub fn required_source_intervals(&self) -> Result<usize, &'static str> {
        if self.rendered >= self.output_samples {
            return Err("bandlimited observation output window is exhausted");
        }
        coverage(self.rendered + 1, self.coverage_width)
    }

    /// Accept one simultaneous original interval-mean frame. A refused frame
    /// changes neither history nor accepted count; channel values retain SI gain.
    pub fn push_frame(&mut self, input: &[f64]) -> Result<(), &'static str> {
        if input.len() != self.report.channels || input.iter().any(|x| !x.is_finite()) {
            return Err("bandlimited observation needs one complete finite source frame");
        }
        if self.accepted >= self.required_source_intervals()? {
            return Err("bandlimited observation source exceeds the next physical output window");
        }
        let begin = (self.accepted % self.report.history_frames) * self.report.channels;
        self.history[begin..begin + self.report.channels].copy_from_slice(input);
        self.accepted += 1;
        Ok(())
    }

    /// Emit one delayed, phase-aligned multichannel frame. Refusal preserves
    /// the caller's output and output cursor. Accepted physical input remains
    /// accepted: this method never rolls back the source that produced it.
    pub fn render_next(&mut self, output: &mut [f64]) -> Result<(), &'static str> {
        if output.len() != self.report.channels {
            return Err("bandlimited observation output has the wrong channel count");
        }
        if self.accepted != self.required_source_intervals()? {
            return Err("bandlimited observation lacks original physical-window coverage");
        }
        let h = self.report.half_support;
        // A source center is i+1/2 and the retained delay is H+1/2 input
        // periods: i_query = output_center/dt - H - 1.
        let coordinate = (self.rendered as f64 + 0.5) / self.physical_width - h as f64 - 1.0;
        let anchor = coordinate.floor();
        let table_phase = (coordinate - anchor) * self.policy.phase_count as f64;
        let phase = (table_phase.floor() as usize).min(self.policy.phase_count - 1);
        let blend = table_phase - phase as f64;
        let start = anchor as i64 - h as i64;
        let end = start + self.report.taps as i64;
        let first = start.max(0) as usize;
        let last = end.max(0) as usize;
        if last > self.accepted
            || (first < last && first < self.accepted.saturating_sub(self.report.history_frames))
        {
            return Err("bandlimited observation stencil is outside accepted retained history");
        }
        self.output.fill(0.0);
        let taps = self.report.taps;
        let channels = self.report.channels;
        let a = &self.coefficients[phase * taps..(phase + 1) * taps];
        let b = &self.coefficients[(phase + 1) * taps..(phase + 2) * taps];
        let mut slot = first % self.report.history_frames;
        for index in first..last {
            let tap = (index as i64 - start) as usize;
            let weight = a[tap] + blend * (b[tap] - a[tap]);
            let frame = &self.history[slot * channels..(slot + 1) * channels];
            for (value, input) in self.output.iter_mut().zip(frame) {
                *value += weight * input;
            }
            slot += 1;
            if slot == self.report.history_frames {
                slot = 0;
            }
        }
        if self.output.iter().any(|x| !x.is_finite()) {
            return Err("bandlimited observation produced a non-finite output frame");
        }
        output.copy_from_slice(&self.output);
        self.rendered += 1;
        Ok(())
    }
}

fn zeros(count: usize) -> Result<Vec<f64>, &'static str> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| "bandlimited observation allocation failed")?;
    values.resize(count, 0.0);
    Ok(values)
}

fn kernel(offset: f64, half: f64, cutoff: f64) -> f64 {
    if offset.abs() >= half {
        return 0.0;
    }
    let phase = PI * offset / half;
    let window = 0.35875
        + 0.48829 * det::cos(phase)
        + 0.14128 * det::cos(2.0 * phase)
        + 0.01168 * det::cos(3.0 * phase)
        - (0.35875 - 0.48829 + 0.14128 - 0.01168);
    let sinc = if offset == 0.0 {
        2.0 * cutoff
    } else {
        det::sin(2.0 * PI * cutoff * offset) / (PI * offset)
    };
    sinc * window
}

// Same endpoint products as the original conservative IntervalClock. Division
// rounding may land on either side of an integer, so verify minimal coverage.
fn coverage(output_end: u64, width: f64) -> Result<usize, &'static str> {
    let end = output_end as f64;
    let count = (end / width).ceil();
    if !count.is_finite() || count < 1.0 || count > EXACT_INTEGER_LIMIT as f64 {
        return Err("bandlimited observation physical interval count is not representable");
    }
    let mut count = count as u64;
    if count > 1 && (count - 1) as f64 * width >= end {
        count -= 1;
    }
    if (count as f64) * width < end {
        count += 1;
    }
    if count > EXACT_INTEGER_LIMIT
        || (count as f64) * width < end
        || (count > 1 && (count - 1) as f64 * width >= end)
    {
        return Err("bandlimited observation cannot distinguish adjacent physical intervals");
    }
    usize::try_from(count)
        .map_err(|_| "bandlimited observation source count exceeds addressable data")
}

#[cfg(test)]
#[path = "resample/tests.rs"]
mod tests;
