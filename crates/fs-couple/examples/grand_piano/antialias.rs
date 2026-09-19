//! Causal integer-rate observation, never a force filter or a material law.
//! The Euler-film resampler requires a complete horizon and reflected edges;
//! this streaming image instead starts from rest and retains explicit latency.
//! Blackman-Harris windowed sinc; power-of-two rates use sparse half-band
//! cascades. Other integer ratios use a direct polyphase output evaluation.
//! Passband is 0..0.45 of OUTPUT rate. For half-bands, rejection is specified
//! at 0.55 output rate and above (foldover into the retained passband), not
//! in the transition around Nyquist. No future samples or delay compensation.
use fs_math::det;
use std::f64::consts::PI;

#[derive(Debug)]
struct Stage {
    ratio: usize, channels: usize, half: usize, center: f64,
    pairs: Vec<(usize, f64)>, history: Vec<f64>, head: usize,
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
            history: vec![0.0; 2 * half * channels], head: 0 }
    }
    // Negative indices address committed history; nonnegative ones address
    // this candidate input block. Preview never writes the delay line.
    fn frame<'a>(&'a self, input: &'a [f64], index: isize) -> &'a [f64] {
        let (data, frame) = if index >= 0 { (input, index as usize) } else {
            let frames = 2 * self.half;
            (&self.history[..], (self.head + frames - index.unsigned_abs()) % frames)
        };
        &data[frame * self.channels..(frame + 1) * self.channels]
    }
    fn preview(&self, input: &[f64], output: &mut [f64]) {
        for (frame, out) in output.chunks_exact_mut(self.channels).enumerate() {
            let end = ((frame + 1) * self.ratio - 1) as isize;
            let center = self.frame(input, end - self.half as isize);
            for (value, x) in out.iter_mut().zip(center) { *value = self.center * x; }
            for &(lag, weight) in &self.pairs {
                let a = self.frame(input, end - lag as isize);
                let b = self.frame(input, end - (2 * self.half - lag) as isize);
                for ((value, a), b) in out.iter_mut().zip(a).zip(b) {
                    *value += weight * (a + b);
                }
            }
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
#[derive(Debug)]
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
        for &(i, w) in &s.pairs { h[i] = w; h[2 * s.half - i] = w; }
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
}
