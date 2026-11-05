//! Cold-planned uniform partitioned overlap-add convolution. Physical travel
//! delay pays for input collection: P <= first tap, so no output shift occurs.
//! Work is bounded by P<=128, <=128 modes and <=4096 retained kernel frames.
//!
//! Inputs and kernels are real, so every spectrum is conjugate-symmetric: only
//! bins 0..=P are kept and multiplied, and two modes share one forward
//! transform (one as the real part, one as the imaginary part), separated by
//! that symmetry. The summed half spectrum is completed before the inverse.
//! Kernels and spectra are stored as separate real and imaginary arrays, so
//! the multiply-accumulate vectorizes without lane shuffles.
// Owner-authorized backport of FrankenPiano receiver_fft (daf1eb708131).
// Transform planning/execution belongs to fs-fft's bounded SmallFft owner.
use fs_fft::{C64 as Complex, SmallFft};
#[derive(Clone)]
pub(crate) struct Uniform {
    plan: SmallFft,
    partition: usize,
    modes: usize,
    parts: usize,
    delay: usize,
    kernel_re: Vec<f64>,
    kernel_im: Vec<f64>,
    spectra_re: Vec<f64>,
    spectra_im: Vec<f64>,
    input: Vec<f64>,
    work: Vec<Complex>,
    acc_re: Vec<f64>,
    acc_im: Vec<f64>,
    sum: Vec<Complex>,
    output: Vec<f64>,
    cursor: usize,
    slot: usize,
    head: usize,
}
impl Uniform {
    pub fn new(rows: &[(usize, Vec<f64>)], modes: usize, partition: usize) -> Result<Self, String> {
        if !(8..=128).contains(&partition)
            || !partition.is_power_of_two()
            || !(1..=128).contains(&modes)
        {
            return Err(
                "FFT receiver requires power-of-two partition 8..128 and 1..128 modes".into(),
            );
        }
        let delay = rows.first().ok_or("empty FFT radiation kernel")?.0;
        let last = rows.last().ok_or("empty FFT radiation kernel")?.0;
        if delay < partition
            || last < delay
            || last - delay >= 4096
            || rows.windows(2).any(|w| w[0].0 >= w[1].0)
            || rows
                .iter()
                .any(|(_, r)| r.len() != modes || r.iter().any(|x| !x.is_finite()))
        {
            return Err("FFT partition exceeds physical delay or kernel budget".into());
        }
        let parts = (last - delay + 1).div_ceil(partition);
        let n = 2 * partition;
        let bins = partition + 1;
        let plan = SmallFft::new(n);
        let mut full = vec![Complex::default(); parts * modes * n];
        for (d, row) in rows {
            let offset = d - delay;
            let part = offset / partition;
            let tap = offset % partition;
            for (mode, &value) in row.iter().enumerate() {
                full[(part * modes + mode) * n + tap].re = value;
            }
        }
        let mut kernel_re = Vec::with_capacity(parts * modes * bins);
        let mut kernel_im = Vec::with_capacity(parts * modes * bins);
        for kernel in full.chunks_exact_mut(n) {
            plan.forward(kernel);
            if kernel
                .iter()
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
            {
                return Err("FFT radiation kernel spectrum overflow".into());
            }
            kernel_re.extend(kernel[..bins].iter().map(|z| z.re));
            kernel_im.extend(kernel[..bins].iter().map(|z| z.im));
        }
        Ok(Self {
            plan,
            partition,
            modes,
            parts,
            delay,
            kernel_re,
            kernel_im,
            spectra_re: vec![0.0; parts * modes * bins],
            spectra_im: vec![0.0; parts * modes * bins],
            input: vec![0.0; partition * modes],
            work: vec![Complex::default(); n],
            acc_re: vec![0.0; bins],
            acc_im: vec![0.0; bins],
            sum: vec![Complex::default(); n],
            output: vec![0.0; delay + n + 1],
            cursor: 0,
            slot: 0,
            head: 0,
        })
    }
    /// Read before committing input: invalid receiver frames must not advance
    /// either convolution history or the following differentiating filter.
    pub fn pressure(&self) -> f64 {
        self.output[self.head]
    }
    pub fn step(&mut self, frame: &[f64]) -> f64 {
        let pressure = self.output[self.head];
        self.output[self.head] = 0.0;
        self.input[self.cursor * self.modes..(self.cursor + 1) * self.modes].copy_from_slice(frame);
        self.cursor += 1;
        if self.cursor == self.partition {
            let n = 2 * self.partition;
            let bins = self.partition + 1;
            for a in (0..self.modes).step_by(2) {
                // Mode a as the real part, mode a + 1 (if any) as the imaginary part.
                let b = (a + 1 < self.modes).then_some(a + 1);
                self.work.fill(Complex::default());
                for (i, z) in self.work[..self.partition].iter_mut().enumerate() {
                    z.re = self.input[i * self.modes + a];
                    if let Some(b) = b {
                        z.im = self.input[i * self.modes + b];
                    }
                }
                self.plan.forward(&mut self.work);
                // X_a[k] = (Z[k] + conj Z[n-k]) / 2 and X_b[k] = (Z[k] - conj Z[n-k]) / 2i.
                let at = (self.slot * self.modes + a) * bins;
                for k in 0..bins {
                    let z = self.work[k];
                    let w = self.work[(n - k) % n];
                    self.spectra_re[at + k] = 0.5 * (z.re + w.re);
                    self.spectra_im[at + k] = 0.5 * (z.im - w.im);
                    if b.is_some() {
                        self.spectra_re[at + bins + k] = 0.5 * (z.im + w.im);
                        self.spectra_im[at + bins + k] = 0.5 * (w.re - z.re);
                    }
                }
            }
            self.acc_re.fill(0.0);
            self.acc_im.fill(0.0);
            let span = self.modes * bins;
            for part in 0..self.parts {
                let previous = (self.slot + self.parts - part) % self.parts;
                let (k, x) = (
                    part * span..(part + 1) * span,
                    previous * span..(previous + 1) * span,
                );
                let h = self.kernel_re[k.clone()]
                    .chunks_exact(bins)
                    .zip(self.kernel_im[k].chunks_exact(bins));
                let x = self.spectra_re[x.clone()]
                    .chunks_exact(bins)
                    .zip(self.spectra_im[x].chunks_exact(bins));
                for ((hr, hi), (xr, xi)) in h.zip(x) {
                    // Complex::mul's terms and the sum's order, per bin.
                    let acc = self.acc_re.iter_mut().zip(self.acc_im.iter_mut());
                    for ((ar, ai), (((&hr, &hi), &xr), &xi)) in
                        acc.zip(hr.iter().zip(hi).zip(xr).zip(xi))
                    {
                        *ar += hr * xr - hi * xi;
                        *ai += hr * xi + hi * xr;
                    }
                }
            }
            for (z, (&re, &im)) in self
                .sum
                .iter_mut()
                .zip(self.acc_re.iter().zip(&self.acc_im))
            {
                *z = Complex { re, im };
            }
            // The output is real: the upper half spectrum mirrors the lower one.
            for k in 1..self.partition {
                let z = self.sum[k];
                self.sum[n - k] = Complex {
                    re: z.re,
                    im: -z.im,
                };
            }
            self.plan.inverse(&mut self.sum);
            let start = self.head + self.delay + 1 - self.partition;
            for (i, z) in self.sum.iter().enumerate() {
                let at = (start + i) % self.output.len();
                self.output[at] += z.re;
            }
            self.cursor = 0;
            self.slot = (self.slot + 1) % self.parts;
        }
        self.head = (self.head + 1) % self.output.len();
        pressure
    }
    #[cfg(test)]
    pub fn state_is_below(&self, limit: f64) -> bool {
        self.input
            .iter()
            .chain(&self.output)
            .chain(&self.spectra_re)
            .chain(&self.spectra_im)
            .chain(&self.acc_re)
            .chain(&self.acc_im)
            .all(|v| v.is_finite() && v.abs() <= limit)
            && self.work.iter().chain(&self.sum).all(|z| {
                z.re.is_finite() && z.im.is_finite() && z.re.abs() <= limit && z.im.abs() <= limit
            })
    }
    #[cfg(test)]
    pub fn flush(&mut self) {
        self.input.fill(0.0);
        self.output.fill(0.0);
        self.spectra_re.fill(0.0);
        self.spectra_im.fill(0.0);
        self.acc_re.fill(0.0);
        self.acc_im.fill(0.0);
        self.work.fill(Complex::default());
        self.sum.fill(Complex::default());
        self.cursor = 0;
        self.slot = 0;
        self.head = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn delayed_partitioned_convolution_matches_direct_at_every_sample_and_wrap() {
        // Odd mode counts leave the last mode without a transform partner.
        for modes in [1, 2, 3] {
            for partition in [8, 16, 32, 64, 128] {
                for delay in [partition, partition + 11] {
                    let rows: Vec<_> = (delay..delay + 3 * partition + 7)
                        .map(|d| {
                            let row = (0..modes)
                                .map(|m| (d as f64 * (m + 1) as f64).sin() * (0.1 + 0.2 * m as f64))
                                .collect();
                            (d, row)
                        })
                        .collect();
                    let mut fft = Uniform::new(&rows, modes, partition).unwrap();
                    let mut history = vec![vec![0.0; modes]; 4096];
                    let (mut error, mut signal) = (0.0, 0.0);
                    for i in 0..4096 {
                        let input: Vec<f64> = (0..modes)
                            .map(|m| (i as f64 * (0.137 + 0.22 * m as f64)).cos())
                            .collect();
                        history[i] = input.clone();
                        let direct: f64 = rows
                            .iter()
                            .filter(|(d, _)| *d <= i)
                            .map(|(d, h)| {
                                h.iter()
                                    .zip(&history[i - d])
                                    .map(|(h, x)| h * x)
                                    .sum::<f64>()
                            })
                            .sum();
                        let candidate = fft.step(&input);
                        error += (direct - candidate).powi(2);
                        signal += direct.powi(2);
                        assert!(
                            (direct - candidate).abs() < 1e-10,
                            "M{modes} P{partition} D{delay} sample{i}"
                        );
                    }
                    assert!(10.0 * (error / signal).log10() < -100.0);
                    fft.flush();
                    assert!(fft.state_is_below(0.0));
                }
            }
        }
        assert!(Uniform::new(&[(127, vec![1.0])], 1, 128).is_err());
    }
}
