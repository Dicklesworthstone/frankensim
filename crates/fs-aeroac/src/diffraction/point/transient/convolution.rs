//! Cached linear convolution of the original real causal basis coefficients.
//! The original FFT owns transforms; this module only pads, normalizes, and
//! contracts the physical Q and dt*Qdot channels. No acoustic frequency fitting.
use super::{AeroacError, Dd, TimedCompactSource, channels, invalid, poll};
use fs_fft::{C64, GeneralFft, GeneralFftBudget};

#[derive(Clone, Debug)]
pub(super) struct Convolution {
    fft: GeneralFft,
    spectra: Vec<Vec<C64>>,
    scales: [f64; 8],
}
impl Convolution {
    pub(super) fn layout(samples: usize) -> Result<(usize, usize), AeroacError> {
        let length = samples
            .checked_mul(2)
            .and_then(|v| v.checked_sub(1))
            .and_then(usize::checked_next_power_of_two)
            .ok_or_else(|| invalid("causal linear convolution size overflow"))?;
        let workspace = GeneralFft::required_workspace_bytes(length).map_err(AeroacError::Fft)?;
        Ok((length, workspace))
    }
    pub(super) fn prepare(
        coefficients: &[[f64; 8]],
        length: usize,
        workspace: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        let fft = GeneralFft::new(
            length,
            GeneralFftBudget {
                max_samples: length,
                max_workspace_bytes: workspace,
            },
        )
        .map_err(AeroacError::Fft)?;
        let mut scratch = filled(fft.scratch_len(), C64::default())?;
        let mut scales = [0.0_f64; 8];
        let mut spectra = Vec::new();
        spectra
            .try_reserve_exact(8)
            .map_err(|_| invalid("causal spectrum-bank allocation failed"))?;
        for component in 0..8 {
            poll(cancelled)?;
            for (i, tap) in coefficients.iter().enumerate() {
                if i % 64 == 0 {
                    poll(cancelled)?;
                }
                scales[component] = scales[component].max(tap[component].abs());
            }
            let mut values = filled(length, C64::default())?;
            if scales[component] > 0.0 {
                for (i, tap) in coefficients.iter().enumerate() {
                    if i % 64 == 0 {
                        poll(cancelled)?;
                    }
                    values[i] = C64::new(tap[component] / scales[component], 0.0);
                }
                fft.forward(&mut values, &mut scratch)
                    .map_err(AeroacError::Fft)?;
            }
            poll(cancelled)?;
            spectra.push(values);
        }
        Ok(Self {
            fft,
            spectra,
            scales,
        })
    }
    pub(super) fn apply(
        &self,
        source: &[TimedCompactSource],
        dt: f64,
        amplitudes: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 4]>, AeroacError> {
        poll(cancelled)?;
        if amplitudes == [0.0; 2] {
            return filled(source.len(), [0.0; 4]);
        }
        let length = self.fft.len();
        let mut source_spectra = [
            filled(length, C64::default())?,
            filled(length, C64::default())?,
        ];
        let mut scratch = filled(self.fft.scratch_len(), C64::default())?;
        let mut work = filled(length, C64::default())?;
        let mut sums = filled(source.len(), [Dd::ZERO; 4])?;
        for channel in 0..2 {
            poll(cancelled)?;
            if amplitudes[channel] == 0.0 {
                continue;
            }
            for (i, sample) in source.iter().enumerate() {
                if i % 64 == 0 {
                    poll(cancelled)?;
                }
                let value = if channel == 0 {
                    sample.strength.volume_flow_m3_s
                } else {
                    dt * sample.strength.volume_flow_rate_m3_s2
                };
                source_spectra[channel][i] = C64::new(value / amplitudes[channel], 0.0);
            }
            self.fft
                .forward(&mut source_spectra[channel], &mut scratch)
                .map_err(AeroacError::Fft)?;
            poll(cancelled)?;
        }
        for component in 0..4 {
            let (a, b) = channels(component);
            for (channel, kernel) in [(0, a), (1, b)] {
                poll(cancelled)?;
                if amplitudes[channel] == 0.0 || self.scales[kernel] == 0.0 {
                    continue;
                }
                for i in 0..length {
                    if i % 64 == 0 {
                        poll(cancelled)?;
                    }
                    let x = source_spectra[channel][i];
                    let y = self.spectra[kernel][i];
                    work[i] = C64::new(
                        (Dd::from_f64(x.re) * Dd::from_f64(y.re)
                            - Dd::from_f64(x.im) * Dd::from_f64(y.im))
                        .to_f64(),
                        (Dd::from_f64(x.re) * Dd::from_f64(y.im)
                            + Dd::from_f64(x.im) * Dd::from_f64(y.re))
                        .to_f64(),
                    );
                }
                self.fft
                    .inverse(&mut work, &mut scratch)
                    .map_err(AeroacError::Fft)?;
                poll(cancelled)?;
                for i in 0..source.len() {
                    if i % 64 == 0 {
                        poll(cancelled)?;
                    }
                    let value = product(work[i].re, amplitudes[channel], self.scales[kernel]);
                    sums[i][component] = sums[i][component] + value;
                }
            }
        }
        let mut result = Vec::new();
        result
            .try_reserve_exact(source.len())
            .map_err(|_| invalid("causal convolution result allocation failed"))?;
        for (i, value) in sums.into_iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let value = value.map(Dd::to_f64);
            if value.iter().any(|v| !v.is_finite()) {
                return Err(invalid("causal physical convolution is unrepresentable"));
            }
            result.push(value);
        }
        poll(cancelled)?;
        Ok(result)
    }
}
fn filled<T: Clone>(count: usize, value: T) -> Result<Vec<T>, AeroacError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| invalid("causal convolution allocation failed"))?;
    result.resize(count, value);
    Ok(result)
}
// Multiply largest and smallest magnitudes first. With three finite factors,
// this avoids an otherwise avoidable intermediate underflow or overflow.
fn product(a: f64, b: f64, c: f64) -> Dd {
    let mut values = [a, b, c];
    values.sort_by(|x, y| x.abs().total_cmp(&y.abs()));
    Dd::from_f64(values[0]) * Dd::from_f64(values[2]) * Dd::from_f64(values[1])
}
