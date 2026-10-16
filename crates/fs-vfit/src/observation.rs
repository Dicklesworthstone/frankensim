//! Causal MIMO observation FIR data and replay. These filters never return a
//! reaction force or assert positive realness: observation is outside the loop.
use fs_math::{c64::C64, det};

/// Total retained coefficients and shared input history, each separately capped.
pub const MAX_OBSERVATION_ENTRIES: usize = 4_194_304;
/// Complete input basis cap; never silently truncate modes/apertures.
pub const MAX_OBSERVATION_INPUTS: usize = 4096;
/// Output components cap.
pub const MAX_OBSERVATION_OUTPUTS: usize = 256;

/// One source-to-receiver-component kernel. Fractional propagation belongs to
/// the taps, not to an unaccounted allpass. The integer prefix is exact.
#[derive(Debug, Clone, PartialEq)]
pub struct FirKernel {
    pub input: usize,
    pub output: usize,
    pub integer_delay_samples: usize,
    pub taps: Vec<f64>,
}
impl FirKernel {
    /// DFT convention exp(-i*omega*t), including the integer propagation delay.
    pub fn transfer_at(&self, omega_rad_per_sample: f64) -> Result<C64, FirError> {
        if !omega_rad_per_sample.is_finite() {
            return Err(FirError::Invalid("nonfinite frequency"));
        }
        let mut h = C64::ZERO;
        for (i, &tap) in self.taps.iter().enumerate() {
            let lag = self
                .integer_delay_samples
                .checked_add(i)
                .ok_or(FirError::Invalid("delay length overflow"))?;
            let phase = -omega_rad_per_sample * lag as f64;
            if !phase.is_finite() {
                return Err(FirError::Invalid("phase overflow"));
            }
            h = h + C64::new(det::cos(phase), det::sin(phase)).scale(tap);
        }
        if !h.re.is_finite() || !h.im.is_finite() {
            return Err(FirError::Invalid("nonfinite FIR transfer"));
        }
        Ok(h)
    }
}

/// Typed cold admission or transactional replay refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FirError {
    Invalid(&'static str),
    Cap { requested: usize, cap: usize },
}
impl core::fmt::Display for FirError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "observation FIR: {self:?}")
    }
}
impl std::error::Error for FirError {}

/// Allocation-free complete-shape admission, also usable by offline producers.
pub fn admit_observation_work(
    inputs: usize,
    outputs: usize,
    coefficients: usize,
    history_entries: usize,
) -> Result<(), FirError> {
    if inputs == 0 || outputs == 0 || coefficients == 0 || history_entries == 0 {
        return Err(FirError::Invalid("empty observation operator"));
    }
    for (requested, cap) in [
        (inputs, MAX_OBSERVATION_INPUTS),
        (outputs, MAX_OBSERVATION_OUTPUTS),
        (coefficients, MAX_OBSERVATION_ENTRIES),
        (history_entries, MAX_OBSERVATION_ENTRIES),
    ] {
        if requested > cap {
            return Err(FirError::Cap { requested, cap });
        }
    }
    Ok(())
}

/// Exact real FIR convolution, zero prehistory, fixed summation order. Input
/// histories are shared across all receiver components. No source is omitted.
#[derive(Clone)]
pub struct ObservationFir {
    inputs: usize,
    outputs: usize,
    kernels: Vec<FirKernel>,
    history: Vec<Vec<f64>>,
    cursor: usize,
    length: usize,
}
impl ObservationFir {
    /// Construct from explicit kernels. Duplicate input/output kernels add;
    /// missing pairs are zero only because the producer explicitly supplied
    /// that sparse operator. A BEM bake supplies every pair and component.
    pub fn new(inputs: usize, outputs: usize, kernels: Vec<FirKernel>) -> Result<Self, FirError> {
        let mut length = 0;
        let mut count = 0usize;
        for k in &kernels {
            if k.input >= inputs
                || k.output >= outputs
                || k.taps.is_empty()
                || k.taps.iter().any(|x| !x.is_finite())
            {
                return Err(FirError::Invalid("invalid kernel index or taps"));
            }
            let end = k
                .integer_delay_samples
                .checked_add(k.taps.len())
                .ok_or(FirError::Invalid("delay length overflow"))?;
            length = length.max(end);
            count = count
                .checked_add(k.taps.len())
                .ok_or(FirError::Invalid("tap count overflow"))?;
        }
        let entries = length
            .checked_mul(inputs)
            .ok_or(FirError::Invalid("history shape overflow"))?;
        admit_observation_work(inputs, outputs, count, entries)?;
        Ok(Self {
            inputs,
            outputs,
            kernels,
            history: vec![vec![0.0; length]; inputs],
            cursor: 0,
            length,
        })
    }
    /// Immutable operator data, suitable for content-addressed artifacts.
    pub fn kernels(&self) -> &[FirKernel] {
        &self.kernels
    }
    /// Transactional one-sample observation. A rejected input/output changes no
    /// history, so retry and replay remain identical after an error.
    pub fn step(&mut self, input: &[f64]) -> Result<Vec<f64>, FirError> {
        if input.len() != self.inputs || input.iter().any(|x| !x.is_finite()) {
            return Err(FirError::Invalid("input shape or nonfinite value"));
        }
        let mut output = vec![0.0; self.outputs];
        for kernel in &self.kernels {
            for (i, &tap) in kernel.taps.iter().enumerate() {
                let lag = i + kernel.integer_delay_samples;
                let value = if lag == 0 {
                    input[kernel.input]
                } else {
                    self.history[kernel.input][(self.cursor + self.length - lag) % self.length]
                };
                output[kernel.output] += tap * value;
            }
        }
        if output.iter().any(|x| !x.is_finite()) {
            return Err(FirError::Invalid("output overflow"));
        }
        for (row, &value) in self.history.iter_mut().zip(input) {
            row[self.cursor] = value;
        }
        self.cursor = (self.cursor + 1) % self.length;
        Ok(output)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_causal_multi_input_replay_is_transactional() {
        let kernels = vec![
            FirKernel {
                input: 0,
                output: 0,
                integer_delay_samples: 2,
                taps: vec![1.0, -0.5],
            },
            FirKernel {
                input: 1,
                output: 0,
                integer_delay_samples: 0,
                taps: vec![2.0],
            },
        ];
        let mut a = ObservationFir::new(2, 1, kernels).unwrap();
        let mut b = a.clone();
        assert!(a.step(&[f64::NAN, 0.0]).is_err());
        for n in 0..8 {
            let x = [
                if n == 0 { 1.0 } else { 0.0 },
                if n == 1 { 3.0 } else { 0.0 },
            ];
            let expected = match n {
                1 => 6.0,
                2 => 1.0,
                3 => -0.5,
                _ => 0.0,
            };
            assert_eq!(a.step(&x).unwrap(), vec![expected]);
            assert_eq!(b.step(&x).unwrap(), vec![expected]);
        }
    }
    #[test]
    fn g0_observation_caps_and_plus_one() {
        assert!(
            admit_observation_work(4096, 256, MAX_OBSERVATION_ENTRIES, MAX_OBSERVATION_ENTRIES)
                .is_ok()
        );
        for request in [
            (4097, 256, 1, 1),
            (4096, 257, 1, 1),
            (1, 1, MAX_OBSERVATION_ENTRIES + 1, 1),
            (1, 1, 1, MAX_OBSERVATION_ENTRIES + 1),
        ] {
            assert!(matches!(
                admit_observation_work(request.0, request.1, request.2, request.3),
                Err(FirError::Cap { .. })
            ));
        }
        assert!(ObservationFir::new(1, 1, vec![]).is_err());
        assert!(
            ObservationFir::new(
                1,
                1,
                vec![FirKernel {
                    input: 1,
                    output: 0,
                    integer_delay_samples: 0,
                    taps: vec![1.0]
                }]
            )
            .is_err()
        );
    }
}
