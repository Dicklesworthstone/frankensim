//! Reusable linear convolution on the original observer's sample history.
//!
//! Overlap-save discards exactly `taps - 1` circularly contaminated samples.
//! The retained block therefore has the same causal indices as scalar FIR
//! stepping, without a block-sized delay or an altered material response.

use super::{AcousticRealizeError, Fft, FftC64, invalid};

/// Maximum prepared observer block, bounding cold FFT and output storage.
pub const MAX_OBSERVER_BLOCK_SAMPLES: usize = 16_384;

#[derive(Debug, Clone)]
pub(super) struct FirBlocks {
    capacity: usize,
    fft: Option<Fft>,
    kernel: Vec<FftC64>,
    spectrum: Vec<FftC64>,
    scratch: Vec<FftC64>,
    output: Vec<f64>,
}

impl FirBlocks {
    pub(super) fn prepare(taps: &[f64], capacity: usize) -> Result<Self, AcousticRealizeError> {
        if capacity == 0 || capacity > MAX_OBSERVER_BLOCK_SAMPLES {
            return Err(invalid(
                "observer block capacity must be in 1..=16384 samples",
            ));
        }
        // The enclosing response owner already admits nonempty finite taps.
        let n = if taps.len() > 32 {
            (taps.len() - 1 + capacity).next_power_of_two()
        } else {
            0
        };
        let fft = (n != 0).then(|| Fft::new(n));
        let mut kernel = vec![FftC64::new(0.0, 0.0); n];
        let mut scratch = kernel.clone();
        if let Some(plan) = &fft {
            for (value, &tap) in kernel.iter_mut().zip(taps) {
                value.re = tap;
            }
            plan.forward(&mut kernel, &mut scratch);
            if kernel
                .iter()
                .any(|value| !value.re.is_finite() || !value.im.is_finite())
            {
                return Err(invalid("observer block kernel is unrepresentable"));
            }
        }
        Ok(Self {
            capacity,
            fft,
            spectrum: vec![FftC64::new(0.0, 0.0); n],
            kernel,
            scratch,
            output: vec![0.0; capacity],
        })
    }

    pub(super) const fn capacity(&self) -> usize {
        self.capacity
    }

    pub(super) fn process(
        &mut self,
        taps: &[f64],
        history: &mut [f64],
        head: &mut usize,
        pressure: &mut [f64],
    ) -> Result<(), AcousticRealizeError> {
        if pressure.len() > self.capacity || pressure.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "observer block needs finite pressure within its prepared capacity",
            ));
        }
        if pressure.is_empty() {
            return Ok(());
        }
        let m = taps.len();
        // Avoid a large FFT for a tiny callback. This is a deterministic work
        // estimate, not a machine-specific benchmark or a real-time guarantee.
        let use_fft = self
            .fft
            .as_ref()
            .is_some_and(|plan| pressure.len() * m > 4 * plan.len() * plan.len().ilog2() as usize);
        if use_fft {
            self.spectrum.fill(FftC64::new(0.0, 0.0));
            // `head` is the next written slot. Skip its oldest sample: only
            // the most recent m-1 inputs precede this block's first sample.
            for (j, value) in self.spectrum[..m - 1].iter_mut().enumerate() {
                value.re = history[(*head + 1 + j) % m];
            }
            for (value, &sample) in self.spectrum[m - 1..].iter_mut().zip(pressure.iter()) {
                value.re = sample;
            }
            let plan = self
                .fft
                .as_ref()
                .expect("FFT path requires a prepared plan");
            plan.forward(&mut self.spectrum, &mut self.scratch);
            for (bin, kernel) in self.spectrum.iter_mut().zip(&self.kernel) {
                *bin = FftC64::new(
                    bin.re * kernel.re - bin.im * kernel.im,
                    bin.re * kernel.im + bin.im * kernel.re,
                );
            }
            plan.inverse(&mut self.spectrum, &mut self.scratch);
            for (out, value) in self
                .output
                .iter_mut()
                .zip(&self.spectrum[m - 1..m - 1 + pressure.len()])
            {
                *out = value.re;
            }
        } else {
            for (i, out) in self.output[..pressure.len()].iter_mut().enumerate() {
                *out = 0.0;
                for (j, &tap) in taps.iter().enumerate() {
                    let sample = if j <= i {
                        pressure[i - j]
                    } else {
                        history[(*head + m - (j - i)) % m]
                    };
                    *out += tap * sample;
                }
            }
        }
        if self.output[..pressure.len()].iter().any(|x| !x.is_finite()) {
            return Err(invalid("observer block output is unrepresentable"));
        }
        // Only accepted input advances the shared scalar/block history. FFT
        // scratch is disposable; a refusal cannot poison the next callback.
        for &sample in pressure.iter() {
            history[*head] = sample;
            *head = if *head + 1 == m { 0 } else { *head + 1 };
        }
        pressure.copy_from_slice(&self.output[..pressure.len()]);
        Ok(())
    }
}

#[cfg(test)]
mod tests;
