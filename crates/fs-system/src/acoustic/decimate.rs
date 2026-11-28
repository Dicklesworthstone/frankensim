//! Causal Blackman-Harris cascade decimator for multirate Rayleigh observation (`.11.9.5`).
//! Matches `fs_couple::pcm_wav::decimate::Decimator::new` (non-differentiating cascade).

use std::f64::consts::PI;

use fs_math::det;

/// Maximum modal channels admitted by the streaming decimator and Rayleigh observer.
pub const MAX_MODES: usize = 128;

#[derive(Clone, Debug)]
struct DecStage {
    ratio: usize,
    ch: usize,
    half: usize,
    center: f64,
    pairs: Vec<(usize, f64)>,
    hist: Vec<f64>,
    head: usize,
}

impl DecStage {
    #[allow(clippy::cast_precision_loss)]
    fn new(ratio: usize, ch: usize, half: usize, hb: bool) -> Self {
        let cut = if hb { 0.25 } else { 0.475 / ratio as f64 };
        let c = |o: usize| {
            let (x, ph) = (o as f64, PI * o as f64 / half as f64);
            let w = 0.35875
                + 0.48829 * det::cos(ph)
                + 0.14128 * det::cos(2.0 * ph)
                + 0.01168 * det::cos(3.0 * ph);
            let s = match o {
                0 => 2.0 * cut,
                _ if hb && o.is_multiple_of(2) => 0.0,
                _ => det::sin(2.0 * PI * cut * x) / (PI * x),
            };
            s * w
        };
        let norm = c(0) + 2.0 * (1..=half).map(c).sum::<f64>();
        let pairs = (0..half)
            .filter_map(|l| Some((l, c(half - l) / norm)).filter(|(_, v)| *v != 0.0))
            .collect();
        let (center, hist, head) = (c(0) / norm, vec![0.0; 2 * half * ch], 0);
        Self {
            ratio,
            ch,
            half,
            center,
            pairs,
            hist,
            head,
        }
    }

    fn frame<'a>(&'a self, inp: &'a [f64], i: isize) -> &'a [f64] {
        let b = i.unsigned_abs();
        let (d, f) = if i >= 0 {
            (inp, b)
        } else {
            (
                &self.hist[..],
                if self.head >= b {
                    self.head - b
                } else {
                    self.head + 2 * self.half - b
                },
            )
        };
        &d[f * self.ch..(f + 1) * self.ch]
    }

    #[allow(clippy::cast_possible_wrap)]
    fn preview(&self, inp: &[f64], out: &mut [f64]) {
        let h = self.half as isize;
        for (fr, o) in out.chunks_exact_mut(self.ch).enumerate() {
            let end = ((fr + 1) * self.ratio - 1) as isize;
            for (v, &x) in o.iter_mut().zip(self.frame(inp, end - h)) {
                *v = self.center * x;
            }
            for &(lag, w) in &self.pairs {
                let l = lag as isize;
                let (a, b) = (self.frame(inp, end - l), self.frame(inp, end - (2 * h - l)));
                for ((v, &av), &bv) in o.iter_mut().zip(a).zip(b) {
                    *v += w * (av + bv);
                }
            }
        }
    }

    fn commit(&mut self, inp: &[f64]) {
        for fr in inp.chunks_exact(self.ch) {
            let off = self.head * self.ch;
            self.hist[off..off + self.ch].copy_from_slice(fr);
            self.head = (self.head + 1) % (2 * self.half);
        }
    }
}

/// Zero-allocation causal Blackman-Harris cascade decimator (`ratio` in `1..=16`).
#[derive(Clone, Debug)]
pub struct Decimator {
    stages: Vec<DecStage>,
    buffers: Vec<Vec<f64>>,
    ratio: usize,
    delay: f64,
}

impl Decimator {
    /// Construct a cold cascade decimator for `ratio` in `1..=16` and `ch` in `1..=128`.
    /// # Errors
    /// Ratio outside `1..=16` or channel count outside `1..=128`.
    #[allow(clippy::cast_precision_loss)]
    pub fn new(ratio: usize, ch: usize) -> Result<Self, &'static str> {
        if !(1..=16).contains(&ratio) || !(1..=MAX_MODES).contains(&ch) {
            return Err("decimator requires ratio in 1..=16 and modes in 1..=128");
        }
        let (mut wk, mut fac) = ([0usize; 17], [0usize; 17]);
        for r in 2..=ratio {
            (wk[r], fac[r]) = (if r == 2 { 41 } else { 80 * r + 1 }, r);
            for f in (2..r).filter(|&f| r.is_multiple_of(f)) {
                let c = (r / f) * (if f == 2 { 9 } else { 8 * f + 1 }) + wk[r / f];
                if c < wk[r] {
                    (wk[r], fac[r]) = (c, f);
                }
            }
        }
        let (mut stages, mut buffers, mut rem, mut delay) =
            (Vec::new(), vec![vec![0.0; ratio * ch]], ratio, 0.0);
        while rem > 1 {
            let (f, hb) = (fac[rem], fac[rem] == 2);
            let half = match (hb, f == rem) {
                (true, true) => 80,
                (true, false) => 16,
                (false, true) => 80 * f,
                _ => 8 * f,
            };
            delay += half as f64 / rem as f64;
            stages.push(DecStage::new(f, ch, half, hb));
            rem /= f;
            buffers.push(vec![0.0; rem * ch]);
        }
        Ok(Self {
            stages,
            buffers,
            ratio,
            delay,
        })
    }

    /// Evaluate candidate decimated frame without mutating delay lines.
    /// # Errors
    /// Non-finite filter output.
    pub fn preview(&mut self, inp: &[f64]) -> Result<(), &'static str> {
        self.buffers[0].copy_from_slice(inp);
        for (i, st) in self.stages.iter().enumerate() {
            let (bef, aft) = self.buffers.split_at_mut(i + 1);
            st.preview(&bef[i], &mut aft[0]);
            if aft[0].iter().any(|v| !v.is_finite()) {
                return Err("decimator output overflow");
            }
        }
        Ok(())
    }

    /// Candidate output slice produced by the last [`Self::preview`].
    #[must_use]
    pub fn pending(&self) -> &[f64] {
        self.buffers.last().map_or(&[], Vec::as_slice)
    }
    /// Commit the last previewed block into the stage delay lines.
    pub fn commit(&mut self) {
        for (s, i) in self.stages.iter_mut().zip(&self.buffers) {
            s.commit(i);
        }
    }
    /// Input mechanics frames per output audio frame (`1..=16`).
    #[must_use]
    pub const fn ratio(&self) -> usize {
        self.ratio
    }
    /// Linear-phase group delay in output frames.
    #[must_use]
    pub const fn delay_output_frames(&self) -> f64 {
        self.delay
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    #[test]
    #[allow(clippy::cast_precision_loss)]
    fn decimator_passband_identity_stopband_rejection_and_refusals() {
        assert!(Decimator::new(0, 1).is_err() && Decimator::new(17, 1).is_err());
        let mut d = Decimator::new(4, 2).unwrap();
        assert!((d.delay_output_frames() - 44.0).abs() < 1e-12 && d.ratio() == 4);
        let (mut pe, mut st, mut blk) = (0.0_f64, 0.0_f64, [0.0; 8]);
        for k in 0..160usize {
            for s in 0..4usize {
                let n = (k * 4 + s) as f64;
                (blk[2 * s], blk[2 * s + 1]) = (
                    det::cos(TAU * 0.1 * n / 4.0),
                    det::cos(TAU * 0.75 * n / 4.0),
                );
            }
            d.preview(&blk).unwrap();
            let (p, q) = (d.pending()[0], d.pending()[1]);
            d.commit();
            if k >= 100 {
                let t = (k * 4 + 3) as f64 / 4.0 - d.delay_output_frames();
                (pe, st) = (pe.max((p - det::cos(TAU * 0.1 * t)).abs()), st.max(q.abs()));
            }
        }
        assert!(pe < 1e-4 && st < 1e-4, "pass={pe} stop={st}");
    }
}
