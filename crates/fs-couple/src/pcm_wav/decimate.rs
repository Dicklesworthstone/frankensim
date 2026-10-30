//! Causal integer-rate observation, never a force filter or a material law.
//! The Euler-film resampler requires a complete horizon and reflected edges;
//! this streaming image instead starts from rest and retains explicit latency.
//! Blackman-Harris windowed sinc; power-of-two rates use sparse half-band
//! cascades. Other integer ratios use a direct polyphase output evaluation.
//! Passband is 0..0.45 of OUTPUT rate. For half-bands, rejection is specified
//! at 0.55 output rate and above (foldover into the retained passband), not
//! in the transition around Nyquist. No future samples or delay compensation.
//! The optional differentiating image observes velocity as acceleration. Its
//! derivative is convolved with the final rate-conversion kernel while cold,
//! not evaluated as a lossy two-sample difference before anti-alias filtering.
use fs_math::det;
use std::f64::consts::PI;

#[derive(Clone, Debug)]
struct Stage {
    ratio: usize, channels: usize, half: usize, center: f64,
    pairs: Vec<(usize, f64)>, history: Vec<f64>, head: usize,
    antisymmetric: bool,
}
impl Stage {
    fn new(ratio: usize, channels: usize, half: usize, halfband: bool) -> Self {
        let cutoff = if halfband { 0.25 } else { 0.475 / ratio as f64 };
        let coefficient = |offset: usize| {
            let x = offset as f64;
            let phase = PI * x / half as f64;
            let window = 0.35875 + 0.48829 * det::cos(phase)
                + 0.14128 * det::cos(2.0 * phase) + 0.01168 * det::cos(3.0 * phase);
            let sinc = if offset == 0 { 2.0 * cutoff }
                else if halfband && offset % 2 == 0 { 0.0 }
                else { det::sin(2.0 * PI * cutoff * x) / (PI * x) };
            sinc * window
        };
        let norm = coefficient(0) + 2.0 * (1..=half).map(coefficient).sum::<f64>();
        let center = coefficient(0) / norm;
        let pairs = (0..half).filter_map(|lag| {
            let value = coefficient(half - lag) / norm;
            (value != 0.0).then_some((lag, value))
        }).collect();
        Self { ratio, channels, half, center, pairs,
            history: vec![0.0; 2 * half * channels], head: 0, antisymmetric: false }
    }

    /// Cold convolution with a bandlimited derivative at THIS stage's input
    /// clock. sinc'(k)=(-1)^k/k for nonzero integer k. Exact antisymmetry
    /// annihilates DC; the first moment gives SI gain for a velocity ramp.
    /// Removing the window's tiny endpoint pedestal makes it continuous at
    /// its finite support. This is differentiation, not acoustic EQ.
    fn differentiate(&mut self, half: usize, input_rate: f64) {
        let mut original = vec![0.0; 2 * self.half + 1];
        original[self.half] = self.center;
        for &(lag, weight) in &self.pairs {
            original[lag] = weight;
            original[2 * self.half - lag] = weight;
        }
        let mut derivative = vec![0.0; 2 * half + 1];
        let mut moment = 0.0;
        for offset in 1..half {
            let x = offset as f64;
            let phase = PI * x / half as f64;
            let window = 0.35875 + 0.48829 * det::cos(phase)
                + 0.14128 * det::cos(2.0 * phase) + 0.01168 * det::cos(3.0 * phase)
                - (0.35875 - 0.48829 + 0.14128 - 0.01168);
            let weight = if offset % 2 == 0 { -window / x } else { window / x };
            derivative[half - offset] = weight;
            derivative[half + offset] = -weight;
            moment += 2.0 * x * weight;
        }
        for value in &mut derivative { *value *= input_rate / moment; }
        let next_half = self.half + half;
        let mut convolution = vec![0.0; 2 * next_half + 1];
        for (i, a) in original.iter().enumerate() {
            for (j, b) in derivative.iter().enumerate() { convolution[i + j] += a * b; }
        }
        self.half = next_half;
        self.center = 0.0;
        // Paired evaluation cancels a constant velocity exactly, rather than
        // subtracting large nearly equal weighted sums at microphone scale.
        self.pairs = (0..next_half).filter_map(|lag| {
            let weight = 0.5 * (convolution[lag] - convolution[2 * next_half - lag]);
            (weight != 0.0).then_some((lag, weight))
        }).collect();
        self.history = vec![0.0; 2 * next_half * self.channels];
        self.head = 0;
        self.antisymmetric = true;
    }

    // Negative indices address committed history; nonnegative ones address
    // this candidate input block. Preview never writes the delay line.
    // Callers reach back 1..=2*half frames and head < 2*half, so one
    // conditional wrap replaces the integer division (two per tap pair).
    fn history_frame(&self, back: usize) -> usize {
        if self.head >= back { self.head - back } else { self.head + 2 * self.half - back }
    }
    fn frame<'a>(&'a self, input: &'a [f64], index: isize) -> &'a [f64] {
        let (data, frame) = if index >= 0 { (input, index as usize) } else {
            (&self.history[..], self.history_frame(index.unsigned_abs()))
        };
        &data[frame * self.channels..(frame + 1) * self.channels]
    }
    fn preview(&self, input: &[f64], output: &mut [f64]) {
        if self.channels == 1 {
            self.preview_one(input, output);
            return;
        }
        for (frame, out) in output.chunks_exact_mut(self.channels).enumerate() {
            let end = ((frame + 1) * self.ratio - 1) as isize;
            let center = self.frame(input, end - self.half as isize);
            for (value, x) in out.iter_mut().zip(center) { *value = self.center * x; }
            if self.antisymmetric {
                for &(lag, weight) in &self.pairs {
                    let a = self.frame(input, end - lag as isize);
                    let b = self.frame(input, end - (2 * self.half - lag) as isize);
                    for ((value, a), b) in out.iter_mut().zip(a).zip(b) {
                        *value += weight * (a - b);
                    }
                }
            } else {
                for &(lag, weight) in &self.pairs {
                    let a = self.frame(input, end - lag as isize);
                    let b = self.frame(input, end - (2 * self.half - lag) as isize);
                    for ((value, a), b) in out.iter_mut().zip(a).zip(b) {
                        *value += weight * (a + b);
                    }
                }
            }
        }
    }
    /// One channel (a summed-pressure derivative): the terms of `preview` in
    /// the same order, read as scalars instead of one-element slices, which
    /// cost more than the arithmetic. Bit-identical.
    fn preview_one(&self, input: &[f64], output: &mut [f64]) {
        let at = |index: isize| if index >= 0 { input[index as usize] }
            else { self.history[self.history_frame(index.unsigned_abs())] };
        for (frame, out) in output.iter_mut().enumerate() {
            let end = ((frame + 1) * self.ratio - 1) as isize;
            let mut value = self.center * at(end - self.half as isize);
            if self.antisymmetric {
                for &(lag, weight) in &self.pairs {
                    value += weight * (at(end - lag as isize) - at(end - (2 * self.half - lag) as isize));
                }
            } else {
                for &(lag, weight) in &self.pairs {
                    value += weight * (at(end - lag as isize) + at(end - (2 * self.half - lag) as isize));
                }
            }
            *out = value;
        }
    }
    fn commit(&mut self, input: &[f64]) {
        for frame in input.chunks_exact(self.channels) {
            let offset = self.head * self.channels;
            self.history[offset..offset + self.channels].copy_from_slice(frame);
            self.head = (self.head + 1) % (2 * self.half);
        }
    }
}

/// One input block is exactly `ratio` interleaved mechanics frames; one output
/// is one interleaved audio frame. All buffers are prepared once. Preview and
/// commit let an acoustic observer reject overflow without losing filter time.
#[derive(Clone, Debug)]
pub struct Decimator {
    stages: Vec<Stage>, buffers: Vec<Vec<f64>>, pending: bool,
    ratio: usize, channels: usize, delay: f64,
}
impl Decimator {
    pub fn new(ratio: usize, channels: usize) -> Result<Self, &'static str> {
        if !(1..=16).contains(&ratio) || !(1..=128).contains(&channels) {
            return Err("decimator needs 1..16 mechanics frames and 1..128 channels");
        }
        let mut stages = Vec::new();
        let mut buffers = vec![vec![0.0; ratio * channels]];
        let mut remaining = ratio;
        let mut delay = 0.0;
        while remaining > 1 {
            let (factor, half, halfband) = if ratio.is_power_of_two() {
                (2, if remaining == 2 { 80 } else { 16 }, true)
            } else { (remaining, 80 * remaining, false) };
            delay += half as f64 / remaining as f64;
            stages.push(Stage::new(factor, channels, half, halfband));
            remaining /= factor;
            buffers.push(vec![0.0; remaining * channels]);
        }
        Ok(Self { stages, buffers, pending: false, ratio, channels, delay })
    }

    /// Observe the time derivative of the supplied samples, in units/second.
    /// `output_rate_hz` names the OUTPUT clock; input frames still arrive at
    /// `ratio * output_rate_hz`. Geometry and physical amplitude are untouched.
    ///
    /// Only the last stage is differentiated, preserving upstream sparse
    /// half-bands and avoiding a separate high-rate FIR per input channel.
    /// Extra latency is 48 output frames at ratio one, otherwise 16 divided
    /// by the final stage's ratio. It is included in `delay_output_frames`.
    /// The legacy `new` constructor remains an undifferentiated observer.
    ///
    /// # Errors
    /// Invalid counts, zero clock or an unrepresentable mechanical clock.
    pub fn new_differentiating(ratio: usize, channels: usize, output_rate_hz: u32)
        -> Result<Self, &'static str> {
        let mut result = Self::new(ratio, channels)?;
        if output_rate_hz == 0 || output_rate_hz.checked_mul(ratio as u32).is_none() {
            return Err("differentiating decimator needs a positive representable sample clock");
        }
        if result.stages.is_empty() {
            // Identity before differentiation; never exposed as a zero-history
            // runnable stage. The differentiator prepares its entire support.
            result.stages.push(Stage { ratio: 1, channels, half: 0, center: 1.0,
                pairs: Vec::new(), history: Vec::new(), head: 0, antisymmetric: false });
            result.buffers.push(vec![0.0; channels]);
        }
        let stage = result.stages.last_mut().expect("derivative stage prepared");
        let half = if stage.ratio == 1 { 48 } else { 16 };
        stage.differentiate(half, f64::from(output_rate_hz) * stage.ratio as f64);
        result.delay += half as f64 / stage.ratio as f64;
        Ok(result)
    }

    pub fn input_frames(&self) -> usize { self.ratio }
    pub fn delay_output_frames(&self) -> f64 { self.delay }
    pub fn multiplies_per_output_frame(&self) -> usize {
        self.stages.iter().zip(&self.buffers[1..]).map(|(s, out)|
            (s.pairs.len() + 1) * out.len()).sum()
    }
    pub fn preview(&mut self, input: &[f64]) -> Result<&[f64], &'static str> {
        self.pending = false;
        if input.len() != self.ratio * self.channels || input.iter().any(|v| !v.is_finite()) {
            return Err("decimator input must contain every finite mechanics frame");
        }
        self.buffers[0].copy_from_slice(input);
        for (i, stage) in self.stages.iter().enumerate() {
            let (before, after) = self.buffers.split_at_mut(i + 1);
            stage.preview(&before[i], &mut after[0]);
            if after[0].iter().any(|v| !v.is_finite()) { return Err("decimator output overflow"); }
        }
        self.pending = true;
        Ok(self.buffers.last().expect("input buffer always present"))
    }
    /// The successful, uncommitted candidate for a shared receiver transaction.
    pub(crate) fn pending_output(&self) -> &[f64] {
        assert!(self.pending, "pending output requires a successful preview");
        self.buffers.last().expect("input buffer always present")
    }
    /// Retained floating-point filter/history values, excluding object metadata.
    pub(crate) fn retained_values(&self) -> usize {
        self.buffers.iter().map(Vec::len).sum::<usize>()
            + self.stages.iter().map(|stage| stage.history.len() + stage.pairs.len()).sum::<usize>()
    }
    pub fn commit(&mut self) {
        assert!(self.pending, "commit requires a successful preview");
        for (stage, input) in self.stages.iter_mut().zip(&self.buffers) { stage.commit(input); }
        self.pending = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn full_taps(s: &Stage) -> Vec<f64> {
        let mut h = vec![0.0; 2 * s.half + 1]; h[s.half] = s.center;
        for &(i, w) in &s.pairs {
            h[i] = w;
            h[2 * s.half - i] = if s.antisymmetric { -w } else { w };
        }
        h
    }
    #[test]
    fn staged_stream_matches_independent_full_convolution_at_every_ratio() {
        for ratio in 1..=16 {
            let channels = 3;
            let input: Vec<f64> = (0..ratio*channels*200).map(|i| (i as f64*0.137).sin()).collect();
            let mut d = Decimator::new(ratio, channels).unwrap();
            let mut reference = input.clone();
            for stage in &d.stages {
                let h = full_taps(stage); let mut next = Vec::new();
                for end in (stage.ratio-1..reference.len()/channels).step_by(stage.ratio) {
                    for c in 0..channels {
                        next.push(h.iter().take(end+1).enumerate().map(|(lag,w)|
                            w*reference[(end-lag)*channels+c]).sum::<f64>());
                    }
                }
                reference = next;
            }
            for (frame, x) in input.chunks_exact(ratio*channels).enumerate() {
                let out = d.preview(x).unwrap();
                for (a,b) in out.iter().zip(&reference[frame*channels..(frame+1)*channels]) {
                    assert!((a-b).abs()<2e-14, "ratio {ratio}: {a} vs {b}");
                }
                d.commit();
            }
        }
    }
    #[test]
    fn one_channel_preview_matches_the_multichannel_path_bit_for_bit() {
        for ratio in 1..=16 {
            for differentiating in [false, true] {
                let make = |channels| if differentiating {
                    Decimator::new_differentiating(ratio, channels, 48_000)
                } else { Decimator::new(ratio, channels) };
                let (mut one, mut two) = (make(1).unwrap(), make(2).unwrap());
                for frame in 0..400 {
                    let x: Vec<f64> = (0..ratio).map(|s| ((frame * ratio + s) as f64 * 0.137).sin()).collect();
                    let pair: Vec<f64> = x.iter().flat_map(|&v| [v, -v]).collect();
                    let (a, b) = (one.preview(&x).unwrap()[0], two.preview(&pair).unwrap()[0]);
                    assert_eq!(a.to_bits(), b.to_bits(), "ratio {ratio}, differentiating {differentiating}, frame {frame}");
                    one.commit();
                    two.commit();
                }
            }
        }
    }
    #[test]
    fn rejection_and_uncommitted_preview_do_not_advance_the_filter() {
        let mut a = Decimator::new(4, 2).unwrap(); let mut b = Decimator::new(4, 2).unwrap();
        let pointers: Vec<_> = a.buffers.iter().map(|v| v.as_ptr()).collect();
        for i in 0..400 {
            a.preview(&[999.0;8]).unwrap(); // abandoned candidate
            assert!(a.preview(&[f64::NAN;8]).is_err());
            assert!(a.preview(&[0.0;7]).is_err());
            let x = [i as f64*0.01;8];
            assert_eq!(a.preview(&x).unwrap(), b.preview(&x).unwrap()); a.commit(); b.commit();
        }
        assert_eq!(pointers, a.buffers.iter().map(|v| v.as_ptr()).collect::<Vec<_>>());
    }
    #[test]
    fn unity_gain_delay_and_bypass_are_explicit() {
        for (ratio, delay) in [(1,0.0),(2,40.0),(4,44.0),(8,46.0),(16,47.0),(3,80.0)] {
            let mut d=Decimator::new(ratio,1).unwrap(); assert_eq!(d.delay_output_frames(),delay);
            for frame in 0..400 {
                let y=d.preview(&vec![1.0;ratio]).unwrap()[0]; d.commit();
                if frame>170 {assert!((y-1.0).abs()<1e-12);}
            }
        }
        assert!(Decimator::new(0,1).is_err()); assert!(Decimator::new(4,129).is_err());
    }
    #[test]
    fn filter_rejects_foldover_into_the_retained_audio_band() {
        for ratio in 2..=16 {
            let d=Decimator::new(ratio,1).unwrap();
            let mut ripple=0.0_f64; let mut alias=0.0_f64;
            for index in 0..=4096 {
                let frequency=ratio as f64*0.5*index as f64/4096.0;
                let mut rate=ratio as f64; let mut gain=1.0;
                for stage in &d.stages {
                    let angle=2.0*PI*frequency/rate;
                    let response=stage.center+2.0*stage.pairs.iter().map(|(lag,w)|
                        w*((stage.half-lag) as f64*angle).cos()).sum::<f64>();
                    gain*=response.abs(); rate/=stage.ratio as f64;
                }
                if frequency<=0.45 {ripple=ripple.max((gain-1.0).abs());}
                if frequency>=0.55 && (frequency-frequency.round()).abs()<=0.45 {
                    alias=alias.max(gain);
                }
            }
            assert!(ripple<0.0001,"ratio {ratio}, ripple {ripple}");
            assert!(alias<0.00001,"ratio {ratio}, foldover {alias}");
        }
    }

    #[test]
    fn differentiating_stream_matches_full_convolution_without_advancing_rejected_trials() {
        for ratio in 1..=16 {
            let channels = 2;
            let input: Vec<f64> = (0..ratio*channels*240)
                .map(|i| (i as f64*0.137).sin()).collect();
            let mut d = Decimator::new_differentiating(ratio, channels, 48_000).unwrap();
            let buffers: Vec<_> = d.buffers.iter().map(|v| v.as_ptr()).collect();
            let histories: Vec<_> = d.stages.iter().map(|s| s.history.as_ptr()).collect();
            let mut reference = input.clone();
            for stage in &d.stages {
                let h = full_taps(stage);
                let mut next = Vec::new();
                for end in (stage.ratio-1..reference.len()/channels).step_by(stage.ratio) {
                    for c in 0..channels {
                        next.push(h.iter().take(end+1).enumerate().map(|(lag,w)|
                            w*reference[(end-lag)*channels+c]).sum::<f64>());
                    }
                }
                reference = next;
            }
            for (frame, x) in input.chunks_exact(ratio*channels).enumerate() {
                d.preview(&vec![999.0; ratio*channels]).unwrap();
                assert!(d.preview(&vec![f64::NAN; ratio*channels]).is_err());
                assert!(d.preview(&x[..x.len()-1]).is_err());
                let out = d.preview(x).unwrap();
                for (a,b) in out.iter().zip(&reference[frame*channels..(frame+1)*channels]) {
                    assert!((a-b).abs()<2e-9, "ratio {ratio}: {a} vs {b}");
                }
                d.commit();
            }
            assert_eq!(buffers, d.buffers.iter().map(|v| v.as_ptr()).collect::<Vec<_>>());
            assert_eq!(histories, d.stages.iter().map(|s| s.history.as_ptr()).collect::<Vec<_>>());
        }
    }

    #[test]
    fn differentiating_gain_is_per_second_and_its_delay_is_explicit() {
        for ratio in 1..=16 {
            for rate in [8_000, 48_000, 192_000] {
                let mut d = Decimator::new_differentiating(ratio, 2, rate).unwrap();
                let base = Decimator::new(ratio, 2).unwrap().delay_output_frames();
                let factor = if ratio.is_power_of_two() && ratio>1 { 2 } else { ratio };
                let extra = if ratio==1 { 48.0 } else { 16.0/factor as f64 };
                assert_eq!(d.delay_output_frames(), base+extra);
                for frame in 0..240 {
                    let mut input = Vec::new();
                    for sub in 0..ratio {
                        let t = (frame*ratio+sub) as f64 / (rate as f64*ratio as f64);
                        input.extend([0.75*t, 1.0]);
                    }
                    let out = d.preview(&input).unwrap();
                    if frame>210 {
                        assert!((out[0]-0.75).abs()<1e-9, "ratio {ratio}, rate {rate}: {}", out[0]);
                        assert_eq!(out[1], 0.0);
                    }
                    d.commit();
                }
            }
        }
        assert!(Decimator::new_differentiating(1,1,0).is_err());
        assert!(Decimator::new_differentiating(2,1,u32::MAX).is_err());
        assert!(Decimator::new_differentiating(17,1,48_000).is_err());
    }

    #[test]
    fn differentiating_passband_and_foldover_keep_the_physical_derivative_gain() {
        for ratio in 1..=16 {
            let d = Decimator::new_differentiating(ratio, 1, 48_000).unwrap();
            let response = |frequency: f64| {
                let mut rate = ratio as f64;
                let mut gain = 1.0;
                for stage in &d.stages {
                    let angle = 2.0*PI*frequency/rate;
                    let value = if stage.antisymmetric {
                        2.0*stage.pairs.iter().map(|(lag,w)|
                            w*((stage.half-lag) as f64*angle).sin()).sum::<f64>()
                    } else {
                        stage.center+2.0*stage.pairs.iter().map(|(lag,w)|
                            w*((stage.half-lag) as f64*angle).cos()).sum::<f64>()
                    };
                    gain *= value;
                    rate /= stage.ratio as f64;
                }
                gain/(2.0*PI*frequency*48_000.0)
            };
            for i in 1..=512 {
                let f = 0.45*i as f64/512.0;
                assert!((response(f)-1.0).abs()<0.0002, "ratio {ratio}, f {f}");
            }
            for i in 1..=1024 {
                let f = 0.5*ratio as f64*i as f64/1024.0;
                if f>=0.55 && (f-f.round()).abs()<=0.45 {
                    assert!(response(f).abs()<0.00002, "ratio {ratio}, alias {f}");
                }
            }
        }
    }
}
